#include <iostream>
#include <vector>
#include <fstream>
#include "../src/include/armadillo"
#include "htslib/sam.h"
#include "../src/saga.h"
#include "../src/interval_tree.h"
#include "../src/gmm.h"
#include "../src/solve_clustering.h"
#include "../src/estimate_error.h"

int main() {
  std::string prefix = "/tmp/var";
  int num_tests = 1;
  int success = 1;

  int32_t primer_offset = 0; 
  uint32_t min_depth = 10;
  uint8_t min_qual = 20;
  uint32_t round_val = 4;
  std::string pair_info = "../data/version_bump_tests/pair_file.tsv";
  std::string bed_file = "../data/version_bump_tests/SARS-CoV-2.primer.bed";

  //TEST 1 - Insertions, deletions, unpaired reads, low quality, paired reads that differ at one site.
  std::string bam_filename = "../data/version_bump_tests/vbump_reads.sorted.bam";
  //std::string var_filename = "../data/version_bump_tests/reads_var.fa";
  std::string ivar_var = "" 
 
  //call variants
  int result = preprocess_reads(bam_filename, bed_file, prefix, "", pair_info, primer_offset, min_depth, min_qual);
  std::vector<variant> base_variants;
  parse_internal_variants(prefix, base_variants, min_depth, round_val, min_qual);
  if(result == 0) success++;
 
  return (num_tests == success) ? 0 : -1;
}
