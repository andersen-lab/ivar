#include "htslib/sam.h"
#include "vector"
#include "string"
#include "ref_seq.h"
#include "interval_tree.h"
#include "allele_functions.h"
#include "site_state.h"
#include "site_aggregator.h"

#ifndef VARIANT_CALLER_H
#define VARIANT_CALLER_H

class variant_caller {
 private:
  static constexpr char seq_nt_lookup[16] = {
      '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
      'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
  };
  uint8_t min_qual;
  uint32_t min_depth;
  ref_antd refantd;
  IntervalTree amplicons;
  site_aggregator sa;
  std::ofstream output_file;
  static const std::string FILE_HEADER;
  static const std::string DELIMITER;

  void get_read_amplicons(uint32_t lower, uint32_t upper, std::vector<ITNode*> &nodes);

 public:
  variant_caller(uint8_t min_qual, uint32_t min_depth, std::string ref_path, std::string gff_path = "");
  bool initialize_region(std::string region);
  ~variant_caller();
  bool parse_read(const bam1_t* read, std::string ref_name, std::vector<site_state> &read_site_states);
  bool parse_paired_reads(const bam1_t* read1, const bam1_t* read2, std::string ref_name, std::vector<site_state> &read_site_states);
  bool parse_single_read(const bam1_t* read, std::string ref_name, std::vector<site_state> &read_site_states) {
    return parse_read(read, ref_name, read_site_states);
  };
  void set_amplicons(IntervalTree &amps);
  void set_refantd(ref_antd &ref);
  void add_variants(std::vector<site_state> &read_site_states);
  void assign_amplicon_to_read(uint32_t lower, uint32_t upper, std::vector<site_state> &read_site_states);
  void merge_reads(std::vector<site_state> &read_site_states_one, std::vector<site_state> &read_site_states_two, std::vector<site_state> &merged_site_states);

  void write_to_file(std::string output_path, std::string ref_name);
};

#endif  // VARIANT_CALLER_H
