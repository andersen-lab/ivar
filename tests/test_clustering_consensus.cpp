#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "htslib/sam.h"
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/ref_seq.h"
#include "../src/parse_gff.h"
#include "../src/call_consensus_clustering.h"
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

  uint32_t min_depth = 5;
  uint8_t min_qual = 20;
  double default_threshold = 0.5;
  double invariant_threshold = 0.99;
  uint32_t n = 2;

  //TEST 1 - manually currated data
  std::string var_filename = "../data/version_bump_tests/vbump_consensus_var.txt";
  std::string consensus_filename = "../data/version_bump_tests/vbump_consensus_ivar.fa";
  std::string reference_file = "../data/version_bump_tests/MN908947.3_sequence.fasta";

  std::vector<double> solution;
  std::vector<double> means;
  std::vector<variant> variants = gmm_model(var_filename, prefix, min_depth, min_qual, solution, means, reference_file, default_threshold, n, invariant_threshold, 1e-3, 1e-2);

  if(!variants.empty()) {
    cluster_consensus(variants, prefix, default_threshold, min_depth, min_qual, solution, means, reference_file);
  }

  std::vector<std::pair<std::string, std::string>> gt_sequences;
  read_consensus(gt_sequences, consensus_filename);
  std::string exp_sequence;
  std::vector<std::pair<std::string, std::string>> exp_sequences;
  read_consensus(exp_sequences, prefix+".fa");

  bool correct = true;
  for (auto itgt = gt_sequences.begin(), itexp = exp_sequences.begin(); itgt != gt_sequences.end() && itexp != exp_sequences.end(); ++itgt, ++itexp) {
    if(itexp->second.size() != itgt->second.size()) {
      correct = false;
      std::cerr << "not same size" << std::endl;
    }
    for(uint32_t i=0; i < itexp->second.size(); i++){
      char a = itgt->second[i];
      char b = itexp->second[i];

      if(a != b){
        std::cerr << "gt " << a << " exp " << b << " " << i+1 << std::endl;
        correct = false;
        break;
      }
    }
  }
  if(correct) success++;

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
