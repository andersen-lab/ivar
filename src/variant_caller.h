#include "htslib/sam.h"
#include "vector"
#include "string"
#include "ref_seq.h"
#include "interval_tree.h"
#include "allele_functions.h"
#include "genomic_position.h"

#ifndef VARIANT_CALLER_H
#define VARIANT_CALLER_H

class variant_caller {
 private:
  static constexpr char seq_nt_lookup[16] = {
      '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
      'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
  };
  uint8_t min_qual;
  ref_antd refantd;
  IntervalTree amplicons;
  std::vector<genomic_position> global_positions;

 public:
  variant_caller(uint8_t min_qual, std::string ref_path, std::string gff_path = "");
  ~variant_caller();
  void parse_read(const bam1_t* read, std::string ref_name, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities);
  void set_amplicons(IntervalTree &amps);
  void add_variants(std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities);
  void get_read_amplicons(uint32_t lower, uint32_t upper, std::vector<ITNode*> &nodes, uint32_t &amp_dist);
  void assign_amplicon_depths(ITNode *node, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities, bool ambiguous);
};

#endif  // IVAR_PARSE_CIGAR_H
