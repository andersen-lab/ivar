#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

int main() {
  int num_tests = 0;
  int success = 0;

  std::string clustering_file = "/tmp/test_default_consensus";


  // 5 positions, each with a single dominant allele: A C G T A.
  // consensus_numbers is set manually here since call_majority_consensus's
  // real callers never populate it (that field is only filled in by
  // assign_variants_solution, which doesn't run on this fallback path).
  std::vector<variant> variants;
  //position 1 has mulitple ambiguity
  {
    variant v{};
    v.position = 1;
    v.nuc = "A";
    v.gapped_freq = 0.80;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  {
    variant v{};
    v.position = 1;
    v.nuc = "C";
    v.gapped_freq = 0.10;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
    {
    variant v{};
    v.position = 1;
    v.nuc = "G";
    v.gapped_freq = 0.10;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  {
    variant v{};
    v.position = 2;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  //position 3 has an ambiguity
  {
    variant v{};
    v.position = 3;
    v.nuc = "A";
    v.gapped_freq = 0.50;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  {
    variant v{};
    v.position = 3;
    v.nuc = "G";
    v.gapped_freq = 0.50;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  {
    variant v{};
    v.position = 4;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  {
    variant v{};
    v.position = 5;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.total_depth = 5;
    v.depth_flag = true;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }

  double default_threshold = 0.0;
  call_majority_consensus(variants, clustering_file, default_threshold);

  // TEST 1 - the written consensus file should contain "AARAN"
  // the first position contains three alleles but the A is the majority, so it is called as A
  // the third position is an ambiguity 0.5/0.5 between AG and should be an R
  // last position below min depth so call an N
  {
    bool pass = true;
    std::ifstream file(clustering_file + "_threshold.fa");
    std::string header, sequence;
    if (std::getline(file, header) && std::getline(file, sequence)) {
      std::string expected = "AARAN";
      if (sequence != expected) {
        pass = false;
        std::cerr << "expected " << expected << " got " << sequence << std::endl;
      }
    } else {
      pass = false;
      std::cerr << "could not read output file " << clustering_file + "_threshold.fa" << std::endl;
    }
    num_tests++;
    if (pass) success++;
  }

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
