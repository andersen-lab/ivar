#include <iostream>
#include <vector>
#include "../src/include/armadillo"
#include "htslib/sam.h"
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/estimate_error.h"
#include "../src/solve_clustering.h"
#include "../src/interval_tree.h"

int main() {
  std::string prefix = "/tmp/consensus";

  //test var
  int num_tests = 2;
  int success = 2;
 
  uint32_t round_val = 4; 
  uint32_t min_depth = 10;
  uint8_t min_qual = 20;
  //std::string var_filename = "../data/version_bump_tests/test_variants.txt";
  std::string var_filename = "../data/version_bump_tests/file_0_var.txt";
   
  std::vector<variant> base_variants; 
  std::vector<variant> variants; 
  uint32_t n = 2;

  //estimation of error
  parse_internal_variants(var_filename, base_variants, min_depth, round_val, min_qual);
  double error_rate = cluster_error(base_variants, min_qual, min_depth);
  float lower_bound = 1-error_rate+0.0001;
  float upper_bound = error_rate-0.0001;
  //parse our test variants file
  uint32_t count = 0;
  set_freq_range_flags(base_variants, lower_bound, upper_bound);
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].outside_freq_range && !base_variants[i].depth_flag){
      variants.push_back(base_variants[i]);
      count++;
    }
  }
  
  //initialize armadillo dataset and populate with frequency data
  arma::mat data(1, count, arma::fill::zeros);

  //(rows, cols) where each columns is a sample
  for(uint32_t i = 0; i < variants.size(); i++){
    double tmp = static_cast<double>(variants[i].gapped_freq);
    data.col(i) = tmp;
  }
  std::vector<double> solution;   
  gaussian_mixture_model retrained = retrain_model(n, data, variants, 2, 0.001);
  solve_clusters(variants, retrained, (double)lower_bound, solution);

  //cluster_consensus(variants, prefix, , double default_t    hreshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means)
  
  exit(0);
  probability_amplicon_frequencies(retrained, base_variants, n);
  
  return (num_tests == success) ? 0 : -1;
}
