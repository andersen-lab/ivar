#include <fstream>
#include <iostream>
#include <string>

#include "../src/allele_functions.h"
#include "../src/call_consensus_pileup.h"

int call_cns_check_outfile(std::string seq_id, std::string prefix, char gap,
                           bool call_min_depth, int min_depth) {
  std::string path = "../data/test.gap.sorted.mpileup";
  std::string expctd_hdr = "";
  std::ifstream mplp(path);
  call_consensus_from_plup(mplp, seq_id, prefix, 20, 0, min_depth, gap,
                           call_min_depth, 1);
  std::ifstream outFile(prefix + ".fa");
  std::string l;
  getline(outFile, l);  // header
  if (seq_id.empty()) {
    char *o = new char[prefix.length() + 1];
    strcpy(o, prefix.c_str());
    expctd_hdr =
        ">Consensus_" + std::string(basename(o)) + "_threshold_0_quality_20";

  } else {
    expctd_hdr = ">" + seq_id;
  }

  return l.compare(expctd_hdr);
}

int main() {
  int num_success = 0;
  std::string c_ =
      "CTGCTGGGTCATGGGCCCATCATGATGGTCTTGGCGATTCTAGCCTTTTTGAGATTCACGGCAATCAAGCCA"
      "TCACTGGGTCTCATCAATAGATGGGGTTCAGTGGGGAAAAAAGAGGCTATGGAAACAATAAAGAAGTTCAAG"
      "AAAGAT------------------------------"
      "AGGAAGGAGAAGAAGAGACGWGGCGCAGATACTAGTGTCGGAATTGTTGGMCTCCTGCTGACCACAGCTATG"
      "GMAGCGGAGGTCACKAGACGTGGGAGTGCATACTATATGTACTTGGACWGAAACGATGCKGGGGAGGCCATA"
      "TCTTTTCCAACCACATTGGGGTTGAATAAGTG";
  num_success =
      call_cns_check_outfile("TESTID", "../data/test.gap", '-', true, 0);
  std::cout << num_success << std::endl;
  num_success += call_cns_check_outfile("", "../data/test.gap", 'N', true, 0);
  std::cout << num_success << std::endl;
  num_success +=
      call_cns_check_outfile("TESTID2", "../data/test.gap", 'N', false, 0);
  std::cout << num_success << std::endl;
  if (num_success == 0) return 0;
  return -1;
}
