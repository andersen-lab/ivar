#include <fstream>
#include <iostream>
#include <string>

#include "../src/allele_functions.h"
#include "../src/call_variants.h"

int call_var_check_outfile(std::string prefix, uint8_t min_qual,
                           uint8_t min_depth, double min_threshold,
                           std::string out[], int len) {
  std::string path = "../data/test.strand.pileup";
  std::ifstream mplp(path);
  call_variants_from_plup(mplp, prefix, min_qual, min_threshold, min_depth,
                          "../data/db/test_ref.fa", "../data/test_strand.gff");
  std::ifstream outFile(prefix + ".tsv");
  std::string l;
  getline(outFile, l);  // Ignore first line
  int comp = 0, ctr = 0;

  while (ctr < len) {
    getline(outFile, l);
    std::cout << l << std::endl;
    std::cout << out[ctr] << " -> CORRECT" << std::endl;
    comp += l.compare(out[ctr]);
    std::cout << l.compare(out[ctr]) << "\n";
    ctr++;
  }
  return comp;
}

int main() {
  int num_success = 0;
  // Quality threshold 20. Frequency threshold: 0.03. Total_DP = 3. Indel passes
  // filters with total_depth 4. Has two lines.
  std::string expected_out[2] = {
      "test\t5\tC\tT\t10\t0\t30\t10\t0\t30\t0.5\t20\t0.000290613\tTRUE\tA:test1\tGCT\tA\tGTT\tV\t2",
      "test\t15\tG\tC\t10\t0\t30\t10\t0\t30\t0.5\t20\t0.000290613\tTRUE\tB:test2\tGTC\tV\tGTG\tV\t2"
      };
  num_success =
      call_var_check_outfile("../data/test_strand", 20, 0, 0.03, expected_out, 2);
  std::cout << num_success << std::endl;

  if (num_success == 0) return 0;
  return -1;
}
