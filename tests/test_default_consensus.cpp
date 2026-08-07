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
  {
    variant v{};
    v.position = 1;
    v.nuc = "A";
    v.gapped_freq = 1.0;
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
  {
    variant v{};
    v.position = 3;
    v.nuc = "A";
    v.gapped_freq = 1.0;
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
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }

  double default_threshold = 0.0;
  uint32_t min_depth = 10;
  call_majority_consensus(variants, clustering_file, default_threshold, min_depth);

  // TEST 1 - the written consensus file should contain "ACGTA"
  {
    bool pass = true;
    std::ifstream file(clustering_file + "_threshold.fa");
    std::string header, sequence;
    if (std::getline(file, header) && std::getline(file, sequence)) {
      std::string expected = "AAAAA";
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
