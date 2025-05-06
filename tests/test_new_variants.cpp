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
#include "../src/allele_functions.h"
#include "../src/parse_gff.h"

int main() {
  std::string prefix = "/tmp/var";
  std::string prefix_2 = "/tmp/var_2";
  int num_tests = 3;
  int success = 0;

  int32_t primer_offset = 0; 
  uint32_t min_depth = 1;
  uint8_t min_qual = 20;
  uint32_t round_val = 4;
  double min_threshold = 0;
  std::string pair_info = "../data/version_bump_tests/pair_file.tsv";
  std::string bed_file = "../data/version_bump_tests/SARS-CoV-2.primer.bed";
  std::string reference_file = "../data/version_bump_tests/MN908947.3_sequence.fasta";

  //TEST 1 - Insertions, deletions, unpaired reads, low quality, paired reads that differ at one site.
  std::string bam_filename = "../data/version_bump_tests/vbump_reads.sorted.bam";
  std::string path = "../data/version_bump_tests/vbump_reads.mpileup";

  //call variants in the new way
  int result = preprocess_reads(bam_filename, bed_file, prefix, "", pair_info, primer_offset, min_depth, min_qual);
  std::vector<variant> new_variants;
  parse_internal_variants(prefix + ".txt", new_variants, min_depth, round_val, min_qual);
  //call variants in the old way
  std::ifstream mplp(path);
  call_variants_from_plup(mplp, prefix, min_qual, min_threshold, min_depth, reference_file, "", true);
  std::vector<variant> old_variants;
  parse_internal_variants(prefix + ".tsv", old_variants, min_depth, round_val, min_qual);

  //compare variants depths between old and new method
  bool depths_match = true;
  for(uint32_t i=0; i < new_variants.size(); i++){
    bool found = false;
    uint32_t position = new_variants[i].position;
    std::string nuc = new_variants[i].nuc; 
    for(uint32_t j=0; j < old_variants.size(); j++){
      if(old_variants[j].position == position && old_variants[j].nuc == nuc){
        found = true;
        if(old_variants[j].depth != new_variants[i].depth){
          std::cerr << "DEPTHS DON'T MATCH POSITION " << position << std::endl;     
          depths_match = false;
          break;
        }
      }
      if(found) break;
    }
  } 
  if(depths_match) success++; 
 
  //compare new variants depths to expected
  bool expected = true;
  for(uint32_t i=0; i < new_variants.size(); i++){
    uint32_t position = new_variants[i].position;
    std::string nuc = new_variants[i].nuc; 
    double total_depth = (double)new_variants[i].total_depth;
    if(position ==58){
      if(total_depth != 6){
        expected = false;
        break;
      }
    } else if(position < 182 || position >= 276 || position == 229){
      if(total_depth != 7){
        expected = false;
        break;
      }
    } else {
      if(total_depth != 8){
        expected = false;
        break;
      }
    }
  } 
  if(expected) success++;

  //TEST 2 - Amplicon flagging correct or incorrect.
  std::string bam_filename_2 = "../data/version_bump_tests/vbump_amplicon.sorted.bam";
  int result_2 = preprocess_reads(bam_filename_2, bed_file, prefix_2, "", pair_info, primer_offset, min_depth, min_qual);
  std::vector<variant> new_variants_2;
  parse_internal_variants(prefix_2 + ".txt", new_variants_2, min_depth, round_val, min_qual);
  bool amp_flags_correct = true;
  for(auto var : new_variants_2){
    if(var.position == 670 && !var.amplicon_flux){
      amp_flags_correct = false;
    }
    if(!var.amplicon_masked){
      amp_flags_correct = false;
    }
  }
  if(amp_flags_correct) success++;  
  std::cerr << "success " << success << " num tests " << num_tests << std::endl; 
  return (num_tests == success) ? 0 : -1;
}
