#include <fstream>
#include <iostream>

#include "../src/get_masked_amplicons.h"

int main() {
  int num_tests = 1, success = 0;
  get_primers_with_mismatches("../data/test.bed", "../data/test.filtered.tsv",
                              "../data/test.masked_primer_indices",
                              "../data/pair_information.tsv");
  std::ifstream masked_indices_file("../data/test.masked_primer_indices.txt");
  std::string indices;
  getline(masked_indices_file, indices);
  if (indices.compare(
          "WNV_400_2_LEFT\tWNV_400_2_RIGHT\tWNV_400_2_LEFT_alt\tWNV_400_1_"
          "LEFT\tWNV_400_1_RIGHT\tWNV_400_1_LEFT_alt\tWNV_400_3_LEFT") == 0)
    success += 1;
  return (num_tests == success) ? 0 : -1;
}
