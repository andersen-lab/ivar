#include "variant_caller.h"

constexpr char variant_caller::seq_nt_lookup[16];

variant_caller::variant_caller(uint8_t min_qual, std::string ref_path, std::string gff_path)
    : min_qual(min_qual),
      refantd(ref_path, gff_path) {}

variant_caller::~variant_caller() {}

void variant_caller::parse_read(const bam1_t* read, std::string ref_name, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities){
  uint32_t total_query_pos=0;
  const uint8_t* query_sequence = bam_get_seq(read);
  std::cout << std::endl;
  uint32_t *cigar = bam_get_cigar(read);
  uint8_t* qual = bam_get_qual(read);
  uint32_t total_ref_pos = read->core.pos + 1; // For 1-based numbering

  for (uint32_t i = 0; i < read->core.n_cigar; i++){
    uint32_t op = bam_cigar_op(cigar[i]);
    uint32_t len = bam_cigar_oplen(cigar[i]);
    if(op == 0){
      for(uint32_t j=0; j < len; j++){
        uint32_t qpos = total_query_pos + j;
        uint32_t rpos = total_ref_pos + j;
        uint8_t tqual = static_cast<uint8_t>(qual[qpos]);
        uint8_t base_code = bam_seqi(query_sequence, qpos);
        char nuc = seq_nt_lookup[base_code];
        positions.push_back(rpos);
        bases.emplace_back(1, nuc);
        qualities.push_back(tqual);
      }
    } else if(op == 2){
      std::string tmp = "-";
      for(uint32_t j=0; j < len; j++){
        tmp += refantd.get_base(total_ref_pos+j, ref_name);
      }
      positions.push_back(total_ref_pos);
      bases.push_back(tmp);
      qualities.push_back(min_qual);
    } else if(op == 1){
      std::string tmp = "+";
      double qual_sum = 0;
      //collect all nucs in insertions
      for(uint32_t j=0; j < len; j++){
        uint32_t qpos = total_query_pos + j;
        uint8_t base_code = bam_seqi(query_sequence, total_query_pos + j);
        char nuc = seq_nt_lookup[base_code];
        tmp += nuc;
        qual_sum += qual[qpos];
      }
      uint8_t avg_qual = static_cast<uint8_t>(qual_sum / len);
      positions.push_back(total_ref_pos-1); // Move insertion position to previous ref
      bases.push_back(tmp);
      qualities.push_back(avg_qual);
    }

    //consumes ref
    if(bam_cigar_type(op) & 2){
      total_ref_pos += len;
    }
    if(bam_cigar_type(op) & 1){
      total_query_pos += len;
    }
  }
}

void variant_caller::set_amplicons(IntervalTree &amps) {
  amplicons = amps;
}

// Unpaired read
void variant_caller::add_variants(std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities){
  if (positions.empty()) return;

  // TODO: Restructure genomic_position so positions with depth = 0, don't need to be stored
  uint32_t max_pos = *std::max_element(positions.begin(), positions.end());
  if (global_positions.size() <= max_pos) {
    global_positions.reserve(max_pos + 1); // reserve prevents reallocations
    while (global_positions.size() <= max_pos) {
      genomic_position tmp;
      tmp.gapped_depth = 0;
      tmp.depth = 0;
      tmp.pos = global_positions.size();
      global_positions.push_back(tmp);
    }
  }

  for (size_t i = 0; i < positions.size(); ++i) {
    const uint32_t pos = positions[i];
    const std::string &base = bases[i];
    const uint32_t qual = qualities[i];

    bool is_del = false;
    bool is_ins = false;

    for (char c : base) {
      if (c == '-') is_del = true;
      else if (c == '+') is_ins = true;
    }

    genomic_position &gpos = global_positions[pos];
    gpos.update_alleles(base, qual);
    if (is_del && !is_ins) {
      //update the gapped depth for all deletion positions
      for(uint32_t k=0; k <bases[i].size()-1; k++){
        global_positions[pos+k].gapped_depth += 1;
      }
    } else if (!is_del && !is_ins) {
      gpos.depth += 1;
      gpos.gapped_depth += 1;
    }
  }
}

void variant_caller::get_read_amplicons(uint32_t start_pos, uint32_t end_pos, std::vector<ITNode*> &nodes, uint32_t &amp_dist){
  amplicons.find_read_amplicon(start_pos, end_pos, nodes);
}

void variant_caller::assign_amplicon_depths(ITNode *node, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities, bool ambiguous = false){
  for(uint32_t i=0; i < positions.size(); i++){
    uint32_t pos = positions[i];
    //for this position, iterate the amplicons associated
    bool found = false;
    bool is_del;
    bool is_ins;
    if(!ambiguous) {
      is_del = bases[i].find('-') != std::string::npos;
      is_ins = bases[i].find('+') != std::string::npos;
    }
    for(uint32_t j=0; j < global_positions[pos].amplicons.size(); j++){
      amplicon_info &amp = global_positions[pos].amplicons[j];
      found = node_compare(node, amp.node);
      if(found){
        if(ambiguous){
          amp.ambiguous_reads += 1;
          break; // TODO: Why break here?
        } else {
          amp.update_alleles(bases[i], qualities[i]);
          if (!is_del && !is_ins) {
            amp.amp_depth += 1;
            //global_positions[pos].gapped_depth += 1;
            //global_positions[pos].depth += 1;
          } else if(!is_ins && is_del){
            amp.amp_depth_gapped += 1;
            for(uint32_t k=0; k < bases[i].size()-1; k++){
              //global_positions[pos+k].gapped_depth += 1;
            }
          }
          break; // TODO: Why break here?
        }
      }
    }
    if(!found){
      //declare a new associated amplicon
      amplicon_info amp;
      amp.node = node;
      amp.amp_depth = 0;
      amp.amp_depth_gapped = 0;
      amp.amp_alleles = populate_basic_alleles();

      if(ambiguous) {
        amp.ambiguous_reads = 1;
      } else {
        if (!is_del && !is_ins) {
          amp.amp_depth += 1;
          //global_positions[pos].gapped_depth += 1;
          //global_positions[pos].depth += 1;
        } else if(!is_ins && is_del){
          amp.amp_depth_gapped += 1;
          for(uint32_t k=0; k <bases[i].size()-1; k++){
            //global_positions[pos+k].gapped_depth += 1;
          }
        }
        amp.update_alleles(bases[i], qualities[i]);
      }
      global_positions[pos].amplicons.push_back(amp);
    }
  }
}