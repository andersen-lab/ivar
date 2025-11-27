#include "variant_caller.h"

constexpr char variant_caller::seq_nt_lookup[16];
const std::string variant_caller::FILE_HEADER="REGION\t"
    "POS\t"
    "REF\t"
    "ALT\t"
    "REF_DP\t"
    "REF_RV\t"
    "REF_QUAL\t"
    "ALT_DP\t"
    "ALT_RV\t"
    "ALT_QUAL\t"
    "ALT_FREQ\t"
    "TOTAL_DP\t"
    "PVAL\t"
    "PASS\t"
    "GFF_FEATURE\t"
    "REF_CODON\t"
    "REF_AA\t"
    "ALT_CODON\t"
    "ALT_AA\t"
    "POS_AA\t"
    "GAPPED_FREQ\t"
    "GAPPED_DEPTH\t"
    "FLAGGED_POS\t"
    "AMP_MASKED\t"
    "STD_DEV\t"
    "AMP_FREQ\t"
    "AMP_NUMBERS\n";
const std::string variant_caller::DELIMITER = "\t";

variant_caller::variant_caller(uint8_t min_qual, std::string ref_path, std::string gff_path)
    : min_qual(min_qual),
      refantd(ref_path, gff_path) {
}

bool variant_caller::initialize_region(std::string region) {
  int64_t ref_len = refantd.get_length(region);
  return sa.initialize(ref_len);
}

variant_caller::~variant_caller() {}

bool variant_caller::parse_read(const bam1_t* read, std::string ref_name, std::vector<site_state> &read_site_states){
  read_site_states.clear();
  read_site_states.reserve(read->core.l_qseq);
  uint32_t total_query_pos=0;
  const uint8_t* query_sequence = bam_get_seq(read);
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
        read_site_states.emplace_back(nuc, tqual, rpos, NUCLEOTIDE);
      }
    } else if(op == 2){
      std::string del_state = "-";
      for(uint32_t j=0; j < len; j++){
        del_state += refantd.get_base(total_ref_pos+j, ref_name);
        if(j > 0){
          // Add gap character for remaining positions
          read_site_states.emplace_back(min_qual, total_ref_pos + j, NUCLEOTIDE);
        }
      }
      read_site_states.emplace_back(del_state, min_qual, total_ref_pos, NUCLEOTIDE);
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
      read_site_states.emplace_back(tmp, avg_qual, total_ref_pos-1, NUCLEOTIDE);
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

bool variant_caller::parse_paired_reads(const bam1_t *read1, const bam1_t *read2, std::string ref_name, std::vector<site_state> &read_site_states) {
  std::vector<site_state> read_site_states_one, read_site_states_two;
  parse_read(read1, ref_name, read_site_states_one);
  parse_read(read2, ref_name, read_site_states_two);
  merge_reads(read_site_states_one, read_site_states_two, read_site_states);
  return true;
}

void variant_caller::set_amplicons(IntervalTree &amps) {
  amplicons = amps;
}

void variant_caller::set_refantd(ref_antd &ref) {
  refantd = ref;
}

void variant_caller::add_variants(std::vector<site_state> &read_site_states){
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
      merged_site_states.emplace_back(read_site_states_one[i]);
      i++;
    } else if (read_site_states_one[i].coordinate > read_site_states_two[j].coordinate) {
      merged_site_states.emplace_back(read_site_states_two[j]);
      j++;
    } else {
      if(read_site_states_one[i].state == read_site_states_two[j].state) {
        read_site_states_one[i].quality = static_cast<uint8_t>(
          static_cast<double>(read_site_states_one[i].quality + read_site_states_two[j].quality) / 2
        );
        merged_site_states.emplace_back(read_site_states_one[i]);
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

void variant_caller::write_to_file(std::string output_path, std::string ref_name) {
  ofstream file;
  file.open(output_path + ".txt", ios::trunc);
  file << FILE_HEADER;
  std::unordered_map<ITNode*, uint32_t> amp_depths;
  std::vector<site_aggregator_stats> aggregated_site_states = sa.get_data();
  // Start iterating from position 1
  for(auto it = aggregated_site_states.begin() + 1; it != aggregated_site_states.end(); ++it) {
    uint32_t coord_position = std::distance(aggregated_site_states.begin(), it);
    const site_coordinate coord = {
      NUCLEOTIDE,
      static_cast<uint32_t>(coord_position)
    };
    const site_aggregator_stats &site_stats = *it;

    // Get ref stats
    const char ref_base = refantd.get_base(coord.position, ref_name);
    const std::string ref_base_str(1, ref_base);
    uint32_t ref_depth = site_stats.get_depth(ref_base_str);
    uint8_t ref_mean_quality = site_stats.get_mean_quality(ref_base_str);

    // Get ref stats for deletion
    const char ref_base_del = refantd.get_base(coord.position - 1, ref_name);
    const std::string ref_base_del_str(1, ref_base_del);
    site_coordinate coord_del;
    coord_del.position = coord.position - 1;
    coord_del.type = coord.type;
    uint32_t ref_depth_del = aggregated_site_states[coord_del.position].get_depth(ref_base_del_str);
    uint32_t ref_mean_quality_del = aggregated_site_states[coord_del.position].get_mean_quality(ref_base_del_str);
    uint32_t total_depth_del = aggregated_site_states[coord_del.position].get_total_depth();
    uint32_t total_gapped_depth_del = aggregated_site_states[coord_del.position].get_total_gapped_depth();

    // Get gapped and ungapped total depth
    uint32_t total_depth = site_stats.get_total_depth();
    uint32_t total_gapped_depth = site_stats.get_total_gapped_depth();
    // Get amplicon depths
    amp_depths.clear();
    if(!sa.calculate_amplicon_depths(coord, amp_depths)){
      std::cerr << "Position " << coord.position << "was not found in " << ref_name << std::endl;
      continue;
    }

    for(auto const &state_stats: site_stats.get_site_state_stats()) {
      const std::string &state = state_stats.get_state();
      if(site_state::is_gap(state))
        continue;

      file << ref_name << DELIMITER; // region
      file << ((site_state::is_deletion(state)) ? coord.position - 1 : coord.position) << DELIMITER; // POS
      file << ((site_state::is_deletion(state)) ? ref_base_del_str : ref_base_str) << DELIMITER;// REF
      file << state << DELIMITER;// ALT
      file << ((site_state::is_deletion(state)) ? ref_depth_del : ref_depth) << DELIMITER; // REF_DP
      file << DELIMITER; // REF_RV
      file << ((site_state::is_deletion(state)) ? static_cast<int>(ref_mean_quality_del) : static_cast<int>(ref_mean_quality)) << DELIMITER; // REF QUAL
      file << state_stats.get_depth() << DELIMITER; // ALT_DP
      file << DELIMITER; // ALT_RV
      file << static_cast<int>(state_stats.get_mean_quality()) << DELIMITER; // ALT_QUAL
      if(site_state::is_deletion(state)){
        file << state_stats.get_depth() / static_cast<double>(total_depth_del) << DELIMITER; // ALT_FREQ
      } else {
        file << state_stats.get_depth() / static_cast<double>(total_depth) << DELIMITER; // ALT_FREQ
      }
      file << total_depth << DELIMITER; // TOTAL_DP
      file << DELIMITER; // PVAL
      file << DELIMITER; // PASS
      file << DELIMITER; // GFF_FEATURE
      file << DELIMITER; //ref codon
      file << DELIMITER; //ref aa
      file << DELIMITER; //alt codon
      file << DELIMITER; //alt aa
      file << DELIMITER; //pos aa
      if(site_state::is_deletion(state)){
        file << state_stats.get_depth() / static_cast<double>(total_gapped_depth_del) << DELIMITER; // GAPPED_FREQ
        file << total_gapped_depth_del << DELIMITER; // GAPPED_DEPTH
      } else {
        file << state_stats.get_depth() / static_cast<double>(total_gapped_depth) << DELIMITER; // GAPPED_FREQ
        file << total_gapped_depth << DELIMITER; // GAPPED_DEPTH
      }

      // Calculating amplicon-specific statistics
      amplicon_variation_data av_data = state_stats.calculate_amplicon_variation(amp_depths);
      file << (av_data.amplicon_frequency_variation ? "True" : "False") << DELIMITER; // FLAGGED_POS -> Variant level
      if(av_data.amplicon_frequency_variation) {
        for (int i = 0; i < av_data.amplicons.size(); i++)
          sa.add_to_masked_amplicons(av_data.amplicons[i]);
      }
      bool is_amplicon_masked = false;
      for(int i = 0; i < av_data.amplicons.size(); i++){
        if (sa.is_amplicon_masked(av_data.amplicons[i])) {
                is_amplicon_masked = true;
                break;
        }
      }
      file << (is_amplicon_masked ? "True" : "False") << DELIMITER; // AMP_MASKED
      file << av_data.stdev << DELIMITER; // STD_DEV
      for(int i = 0; i < av_data.frequencies.size(); i++){
        if(i > 0)
          file << ",";
        file << av_data.frequencies[i]; // AMP FREQ
      }
      file << DELIMITER;
      for(int i = 0; i < av_data.amplicons.size(); i++){
        if(i > 0)
          file << ",";
        file << av_data.amplicons[i]->to_string(); // AMP NUMBERS
      }
      file << "\n";
    }
  }
  file.close();
}