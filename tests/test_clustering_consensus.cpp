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
  //test var
  int num_tests = 2;
  int success = 2;
 
  uint32_t round_val = 4; 
  uint32_t min_depth = 10;
  uint8_t min_qual = 20;
  float lower_bound = 0.03;
  float upper_bound = 0.97;
  std::string var_filename = "../data/version_bump_tests/test_variants.txt";
    
  std::vector<variant> base_variants; 
  uint32_t n = 4;

  //parse our test variants file
  parse_internal_variants(var_filename, base_variants, min_depth, lower_bound, upper_bound, round_val, min_qual);

  //initialize armadillo dataset and populate with frequency data
  arma::mat data(1, base_variants.size(), arma::fill::zeros);

  //(rows, cols) where each columns is a sample
  for(uint32_t i = 0; i < base_variants.size(); i++){
    double tmp = static_cast<double>(base_variants[i].gapped_freq);
    data.col(i) = tmp;
  }

  gaussian_mixture_model retrained = retrain_model(n, data, base_variants, 2, 0.001);
  solve_clusters(base_variants, retrained);
  exit(0);
  for(auto m : retrained.means){
    std::cerr << m << std::endl;
  }

  probability_amplicon_frequencies(retrained, base_variants, n);
  
  return (num_tests == success) ? 0 : -1;
}
