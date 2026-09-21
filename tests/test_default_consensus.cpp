#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

std::string read_sequence(const std::string &path){
  std::ifstream file(path);
  std::string header, sequence;
  if(std::getline(file, header) && std::getline(file, sequence)){
    return sequence;
  }
  return "";
}

void check(const std::string &rule, const std::string &expected, const std::string &actual, int &num_tests, int &success){
  num_tests++;
  if(actual == expected){
    success++;
  } else {
    std::cerr << rule << " - expected " << expected << " got " << actual << std::endl;
  }
}

variant make_variant(uint32_t position, const std::string &nuc, double gapped_freq){
  variant v{};
  v.position = position;
  v.nuc = nuc;
  v.gapped_freq = gapped_freq;
  v.total_depth = 100;
  return v;
}

int main() {
  int num_tests = 0;
  int success = 0;

  std::string clustering_file = "/tmp/test_default_consensus";


  // 5 positions, each with a single dominant allele: A C G T A.
  // consensus_numbers is vestigial here - call_majority_consensus votes over every
  // record regardless of it. The later tests leave it empty on purpose.
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
  //position 2 should still be called as A since it is the majority allele
  {
    variant v{};
    v.position = 2;
    v.nuc = "A";
    v.gapped_freq = 0.95;
    v.total_depth = 100;
    v.consensus_numbers = {0};
    variants.push_back(v);
  }
  {
    variant v{};
    v.position = 2;
    v.nuc = "G";
    v.gapped_freq = 0.05;
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

  default_threshold = 0.90;
  call_majority_consensus(variants, clustering_file, default_threshold);
  // TEST 2 - the written consensus file should contain "VARAN"
  // the first position contains three alleles but all are below the threshold of 0.9, so it is called as an ambiguity V
  // the third position is an ambiguity 0.5/0.5 between AG and should be an R
  // last position below min depth so call an N
  {
    bool pass = true;
    std::ifstream file(clustering_file + "_threshold.fa");
    std::string header, sequence;
    if (std::getline(file, header) && std::getline(file, sequence)) {
      std::string expected = "VARAN";
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

  // Every case below leaves consensus_numbers empty, which is what the real callers
  // produce for all but the invariant records. Deletions are hand built, so position is
  // the first deleted base with no +1 adjustment. write_consensus_to_file strips "-", so
  // a deleted position is a missing character in the expected string.
  double t = 0.75;

  // TEST 3 - records with no consensus_numbers still reach the consensus
  {
    std::vector<variant> v = {make_variant(1, "A", 0.97), make_variant(2, "C", 0.95)};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 3 unassigned records reach the consensus", "AC",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  // TEST 4 - insertion at or above the threshold is written
  {
    std::vector<variant> v = {make_variant(1, "A", 0.95), make_variant(1, "+GAG", 0.90)};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 4 insertion above threshold written", "AGAG",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  // TEST 5 - insertion below the threshold is dropped
  {
    std::vector<variant> v = {make_variant(1, "C", 0.99), make_variant(1, "+TT", 0.20)};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 5 insertion below threshold dropped", "C",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  // TEST 6 - the best supported insertion wins, not the last one read
  {
    std::vector<variant> v = {make_variant(1, "G", 0.99), make_variant(1, "+C", 0.85),
                              make_variant(1, "+AAAA", 0.80), make_variant(1, "+TT", 0.78)};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 6 best supported insertion wins", "GC",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  // TEST 7 - insertions honor the depth and quality filters
  {
    variant blocked_depth = make_variant(1, "+AA", 0.95);
    blocked_depth.depth_flag = true;
    variant blocked_qual = make_variant(2, "+CC", 0.95);
    blocked_qual.qual_flag = true;
    std::vector<variant> v = {make_variant(1, "T", 0.99), blocked_depth,
                              make_variant(2, "A", 0.99), blocked_qual};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 7 insertions honor depth and quality filters", "TA",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  // TEST 8 - a sub threshold deletion must not blank the bases across its span
  {
    std::vector<variant> v = {make_variant(1, "A", 0.97), make_variant(1, "-AA", 0.05),
                              make_variant(2, "C", 0.95)};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 8 sub threshold deletion does not blank its span", "AC",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  // TEST 9 - a dominant deletion still deletes its span
  {
    std::vector<variant> v = {make_variant(1, "G", 0.20), make_variant(1, "-GG", 0.80),
                              make_variant(2, "T", 0.15), make_variant(3, "A", 0.99)};
    call_majority_consensus(v, clustering_file, t);
    check("TEST 9 dominant deletion deletes its span", "A",
          read_sequence(clustering_file + "_threshold.fa"), num_tests, success);
  }

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
