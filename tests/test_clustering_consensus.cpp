#include <iostream>
#include <vector>
#include <fstream>
#include "../src/include/armadillo"
#include "htslib/sam.h"
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/ref_seq.h"
#include "../src/parse_gff.h"
#include "../src/call_consensus_clustering.h"
#include "../src/estimate_error.h"
#include "../src/solve_clustering.h"
#include "../src/interval_tree.h"

void read_consensus(std::vector<std::pair<std::string, std::string>> &sequences, std::string filename){
  std::ifstream file(filename);
  std::string sequence;
  std::string name;
  std::string tmp;
  while (std::getline(file, tmp)) {
    if (tmp.find(">") != std::string::npos){
      if(sequence.size() > 0){
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),[](unsigned char c) { return std::toupper(c); });
        sequences.emplace_back(name, sequence);
      }
      name = tmp;
      sequence.clear();
      continue;
    }
    sequence += tmp;
  }
  if(sequence.size() > 0){
    sequences.emplace_back(name, sequence);
  }  
}

int main() {
  std::string prefix = "/tmp/consensus";
  int num_tests = 1;
  int success = 0;
 
  uint32_t round_val = 4; 
  uint32_t min_depth = 10;
  uint8_t min_qual = 20;
  double default_threshold = 0.5;

  //TEST 1 - manually currated data
  std::string var_filename = "../data/version_bump_tests/test_variants_small.tsv";
  std::string consensus_filename = "../data/version_bump_tests/test_consensus.fa";
 
  std::vector<variant> base_variants; 
  std::vector<variant> variants; 
  uint32_t n = 2;

  //estimation of error
  parse_internal_variants(var_filename, base_variants, min_depth, round_val, min_qual);
  double error_rate = cluster_error(base_variants, min_qual, min_depth);
  double lower_bound = 1-error_rate+0.0001;
  double upper_bound = error_rate-0.0001;
  //parse our test variants file
  uint32_t count = 0;
  set_freq_range_flags(base_variants, lower_bound, upper_bound);
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].outside_freq_range && !base_variants[i].depth_flag){
      variants.push_back(base_variants[i]);
      count++;
    }
  }
  arma::mat data(1, count, arma::fill::zeros);

  //(rows, cols) where each columns is a sample
  for(uint32_t i = 0; i < variants.size(); i++){
    double tmp = static_cast<double>(variants[i].gapped_freq);
    data.col(i) = tmp;
  }
  std::vector<double> solution;   
  gaussian_mixture_model retrained = retrain_model(n, data, variants, 2, 0.001);
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].outside_freq_range || base_variants[i].depth_flag){
      variants.push_back(base_variants[i]);   
    }
  }
  solve_clusters(variants, retrained, lower_bound, solution);
  cluster_consensus(variants, prefix, default_threshold, min_depth, min_qual, solution, retrained.means); 

  std::vector<pair<std::string, std::string>> gt_sequences;
  read_consensus(gt_sequences, consensus_filename);
  std::string exp_sequence;
  std::vector<pair<std::string, std::string>> exp_sequences;
  read_consensus(exp_sequences, prefix+".fa");
  bool correct = true;
  for (auto itgt = gt_sequences.begin(), itexp = exp_sequences.begin(); itgt != gt_sequences.end() && itexp != exp_sequences.end(); ++itgt, ++itexp) {
    for(uint32_t i=0; i < itexp->second.size(); i++){
      char a = itgt->second[i];
      char b = itexp->second[i];
      if(a != b){
        correct = false;
        break;
      }
    }
  }
  
  if(correct) success++;  
  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
