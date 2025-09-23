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
#include "../src/call_variants.h"
#include "../src/ref_seq.h"
#include "../src/call_consensus_clustering.h"
#include "../src/allele_functions.h"
#include "../src/parse_gff.h"
#include "../src/genomic_position.h"

variant find_identical_variant(variant search_var, std::vector<variant> all_variants){
  variant null;
  null.nuc = "None";
  for (auto var : all_variants){
    if(var.position == search_var.position && var.nuc == search_var.nuc) return(var);
  }
  return(null);
}

bool test_depth(std::vector<variant> new_variants, std::vector<variant> old_variants){
  bool depths_match = true;
  for(uint32_t i=0; i < new_variants.size(); i++){
    variant match = find_identical_variant(new_variants[i], old_variants);
    std::cerr << match.nuc << " " << new_variants[i].position << std::endl;
    exit(0);
  }
  return(depths_match);
}

int main() {
  std::string prefix = "/tmp/var";
  int num_tests = 1;
  int success = 0;

  int32_t primer_offset = 0;
  uint32_t min_depth = 1;
  uint8_t min_qual = 20;
  uint32_t round_val = 4;
  double min_threshold = 0;

  std::string pair_info = "../data/version_bump_tests/pair_file.tsv";
  std::string bed_file = "../data/version_bump_tests/SARS-CoV-2.primer.bed";
  std::string reference_file = "../data/version_bump_tests/MN908947.3_sequence.fasta";

  std::string bam_filename = "../data/version_bump_tests/test_amplicon_ambig.bam";

  //call variants
  int result = preprocess_reads(bam_filename, bed_file, prefix, "", pair_info, primer_offset, min_depth, min_qual, reference_file);
  std::vector<variant> new_variants;
  parse_internal_variants(prefix + ".txt", new_variants, min_depth, round_val, min_qual);

  //make sure flags are properly set


  std::cerr << "success " << success << " num tests " << num_tests << std::endl;
  return (num_tests == success) ? 0 : -1;
}
