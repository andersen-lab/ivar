#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <limits>
#include <set>
#include <cstdio>
#include "../src/gmm.h"
#include "../src/call_consensus_clustering.h"

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
  std::string prefix = "/tmp/consensus_integration";
  int num_tests = 1;
  int success = 0;

  uint32_t min_depth = 10;
  uint8_t min_qual = 20;
  double default_threshold = 0.5;
  double invariant_threshold = 0.97;

  uint32_t n = 6;

  //TEST 1 - full consensus against manually currated ground truth
  std::string var_filename = "../../data/version_bump_tests/vbump_mixture_var.txt";
  std::string consensus_filename = "../../data/version_bump_tests/vbump_mixture_ivar.fa";

  double covariance_prior = 1e-3;
  double mean_precision_prior = 1e-2;
  double half_normal_covariance_prior = 1e-3;
  double min_cluster_fraction = 0.05;
  uint32_t min_cluster_points = 0;
  double amplicon_stdev = DEFAULT_AMPLICON_STDEV;
  bool dump_variants = false;
  double freq_outlier_threshold = DEFAULT_FREQ_OUTLIER_THRESHOLD;

  //clear output from any earlier run, so only this run's files can be compared
  std::remove((prefix + ".fa").c_str());
  std::remove((prefix + "_threshold.fa").c_str());
  std::remove((prefix + "_gmm_1d_results.txt").c_str());

  std::vector<double> solution;
  std::vector<double> means;
  std::vector<variant> variants = gmm_model(var_filename, prefix, min_depth, min_qual, solution, means, default_threshold,
                                            n, invariant_threshold, covariance_prior, mean_precision_prior,
                                            half_normal_covariance_prior, min_cluster_fraction, min_cluster_points,
                                            amplicon_stdev, dump_variants, freq_outlier_threshold);
  std::cerr << "effective means:";
  for(double m : means){
    std::cerr << " " << m;
  }
  std::cerr << std::endl;
  if(!variants.empty()) {
    cluster_consensus(variants, prefix, solution);
  }
  std::vector<std::pair<std::string, std::string>> gt_sequences;
  read_consensus(gt_sequences, consensus_filename);
  std::vector<std::pair<std::string, std::string>> exp_sequences;
  read_consensus(exp_sequences, prefix+".fa");

  bool correct = true;
  if(variants.empty()){
    std::cerr << "clustering found no single solution, only the threshold consensus was written" << std::endl;
    correct = false;
  }
  if(gt_sequences.empty() || exp_sequences.empty()){
    std::cerr << "nothing to compare: expected " << gt_sequences.size()
              << " sequences, produced " << exp_sequences.size() << std::endl;
    correct = false;
  } else {
    //compare the expected consensus to the most abundant genome, which is written first. The
    //pileup consensus runs the full reference while the cluster consensus stops at the last
    //variant, so trailing N (no coverage) is ignored on both
    std::string gt = gt_sequences[0].second;
    std::string exp = exp_sequences[0].second;
    gt.erase(gt.find_last_not_of('N') + 1);
    exp.erase(exp.find_last_not_of('N') + 1);
    //genome position 990 is the designed amplicon failure: the pileup consensus calls it M, the
    //cluster consensus N. Consensus position 987, after the 3 bp deletion at 600-602
    std::set<size_t> excluded = {987};
    if(gt.size() != exp.size()){
      correct = false;
      std::cerr << "length differs: expected " << gt.size() << " produced " << exp.size() << std::endl;
    }
    for(size_t i = 0; i < std::min(gt.size(), exp.size()); i++){
      if(gt[i] != exp[i] && excluded.count(i+1) == 0){
        correct = false;
        std::cerr << "consensus position " << i+1 << ": expected " << gt[i] << " produced " << exp[i] << std::endl;
      }
    }
  }
  if(correct) success++;

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
