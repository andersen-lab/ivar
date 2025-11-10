#include "variant_caller.h"

constexpr char variant_caller::seq_nt_lookup[16];

variant_caller::variant_caller(uint8_t min_qual, std::string ref_path, std::string gff_path)
    : min_qual(min_qual),
      refantd(ref_path, gff_path) {}

variant_caller::~variant_caller() {}

bool variant_caller::parse_read(const bam1_t* read, std::string ref_name, std::vector<site_state> &read_site_states){
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
        site_state ss;
        ss.set_nucleotide(std::string(1, nuc), tqual, rpos);
        read_site_states.emplace_back(ss);
      }
    } else if(op == 2){
      std::string tmp = "-";
      for(uint32_t j=0; j < len; j++){
        tmp += refantd.get_base(total_ref_pos+j, ref_name);
        if(j > 0){
          // Add gap character for remaining positions
          site_state ss;
          ss.set_nucleotide_gap(min_qual, total_ref_pos + j);
          read_site_states.emplace_back(ss);
        }
      }
      site_state ss;
      ss.set_nucleotide(tmp, min_qual, total_ref_pos);
      read_site_states.emplace_back(ss);
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
      site_state ss;
      ss.set_nucleotide(tmp, avg_qual, total_ref_pos-1); // Move insertion position to previous ref
      read_site_states.emplace_back(ss);
    }

    //consumes ref
    if(bam_cigar_type(op) & 2){
      total_ref_pos += len;
    }
    if(bam_cigar_type(op) & 1){
      total_query_pos += len;
    }
  }
  return true;
}

void variant_caller::set_amplicons(IntervalTree &amps) {
  amplicons = amps;
}

void variant_caller::set_refantd(ref_antd &ref) {
  refantd = ref;
}

void variant_caller::add_variants(std::vector<site_state> read_site_states){
  sa.aggregate(read_site_states);
}

void variant_caller::get_read_amplicons(uint32_t start_pos, uint32_t end_pos, std::vector<ITNode*> &nodes){
  amplicons.find_read_amplicon(start_pos, end_pos, nodes);
}

void variant_caller::merge_reads(std::vector<site_state> &read_site_states_one, std::vector<site_state> &read_site_states_two, std::vector<site_state> &merged_site_states) {
  size_t i = 0, j = 0;
  merged_site_states.reserve(read_site_states_one.size() + read_site_states_two.size());
  while (i < read_site_states_one.size() && j < read_site_states_two.size()) {
    if (read_site_states_one[i].coordinate < read_site_states_two[j].coordinate) {
      merged_site_states.push_back(read_site_states_one[i]);
      i++;
    } else if (read_site_states_one[i].coordinate > read_site_states_two[j].coordinate) {
      merged_site_states.push_back(read_site_states_two[j]);
      j++;
    } else {
      if(read_site_states_one[i].state == read_site_states_two[j].state) {
        read_site_states_one[i].quality = static_cast<uint8_t>(
          static_cast<double>(read_site_states_one[i].quality + read_site_states_two[j].quality) / 2
        );
        merged_site_states.push_back(read_site_states_one[i]);
      }
      i++;
      j++;
    }
  }

  merged_site_states.insert(merged_site_states.end(), read_site_states_one.begin() + i, read_site_states_one.end());
  merged_site_states.insert(merged_site_states.end(), read_site_states_two.begin() + j, read_site_states_two.end());
}

void variant_caller::assign_amplicon_to_read(uint32_t lower, uint32_t upper, std::vector<site_state> &read_site_states) {
  std::vector<ITNode*> assigned_amplicons;
  get_read_amplicons(lower, upper, assigned_amplicons);
  if (assigned_amplicons.empty() || assigned_amplicons.size() > 1) {
    for(int i = 0; i < read_site_states.size(); i++){
      read_site_states[i].set_amplicon(nullptr, true);
    }
  } else if (assigned_amplicons.size() == 1){
    for(int i = 0; i < read_site_states.size(); i++){
      read_site_states[i].set_amplicon(assigned_amplicons[0], false);
    }
  }
}