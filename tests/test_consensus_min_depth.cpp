#include <fstream>
#include <iostream>
#include <string>

#include "../src/allele_functions.h"
#include "../src/call_consensus_pileup.h"

int call_cns_check_outfile(std::string input_id, std::string prefix,
                           std::string cns, char gap, bool call_min_depth,
                           int min_depth) {
  std::string path = "../data/test.gap.sorted.mpileup";
  std::ifstream mplp(path);
  call_consensus_from_plup(mplp, input_id, prefix, 20, 0, min_depth, gap,
                           call_min_depth, 1);
  std::ifstream outFile(prefix + ".fa");
  std::string l;
  getline(outFile, l);  // Ignore first line
  getline(outFile, l);
  return l.compare(cns);
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
  std::string cN =
      "CTGCTGGGTCATGGGCCCATCATGATGGTCTTGGCGATTCTAGCCTTTTTGAGATTCACGGCAATCAAGCCA"
      "TCACTGGGTCTCATCAATAGATGGGGTTCAGTGGGGAAAAAAGAGGCTATGGAAACAATAAAGAAGTTCAAG"
      "AAAGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGAAGGAGAAGAAGAGACGWGGCGCAGATACTAGT"
      "GTCGGAATTGTTGGMCTCCTGCTGACCACAGCTATGGMAGCGGAGGTCACKAGACGTGGGAGTGCATACTAT"
      "ATGTACTTGGACWGAAACGATGCKGGGGAGGCCATATCTTTTCCAACCACATTGGGGTTGAATAAGTG";
  std::string ck =
      "CTGCTGGGTCATGGGCCCATCATGATGGTCTTGGCGATTCTAGCCTTTTTGAGATTCACGGCAATCAAGCCA"
      "TCACTGGGTCTCATCAATAGATGGGGTTCAGTGGGGAAAAAAGAGGCTATGGAAACAATAAAGAAGTTCAAG"
      "AAAGATAGGAAGGAGAAGAAGAGACGWGGCGCAGATACTAGTGTCGGAATTGTTGGMCTCCTGCTGACCACA"
      "GCTATGGMAGCGGAGGTCACKAGACGTGGGAGTGCATACTATATGTACTTGGACWGAAACGATGCKGGGGAG"
      "GCCATATCTTTTCCAACCACATTGGGGTTGAATAAGTG";
  num_success =
      call_cns_check_outfile("", "../data/test.gap", c_, '-', true, 0);
  std::cout << num_success << std::endl;
  num_success +=
      call_cns_check_outfile("TESTID", "../data/test.gap", cN, 'N', true, 0);
  std::cout << num_success << std::endl;
  num_success +=
      call_cns_check_outfile("TESTID", "../data/test.gap", ck, 'N', false, 0);
  std::cout << num_success << std::endl;
  if (num_success == 0) return 0;
  return -1;
}
