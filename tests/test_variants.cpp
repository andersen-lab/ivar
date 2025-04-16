#include <fstream>
#include <iostream>
#include <string>

#include "../src/allele_functions.h"
#include "../src/call_variants.h"

int call_var_check_outfile(std::string prefix, uint8_t min_qual,
                           uint8_t min_depth, double min_threshold,
                           std::string out[], int len) {
  std::string path = "../../data/test.indel.mpileup";
  std::ifstream mplp(path);
  call_variants_from_plup(mplp, prefix, min_qual, min_threshold, min_depth,
                          "../data/db/test_ref.fa", "../data/test.gff", false);
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
  std::string t_20_02_1[4] = {
      "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0."
      "2\tFALSE\tgeneKatyPerry:id-test3\tAAG\tK\tATG\tM\t70",
      "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0.2\tFALSE\tid-"
      "testedit1\tGAA\tE\tGAT\tD\t70",
      "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0."
      "2\tFALSE\tgeneTaylorSwift:id-testedit2\tAGA\tR\tTGA\t*\t70",
      "test\t210\tA\t+GT\t1\t1\t41\t1\t0\t20\t0.25\t4\t0."
      "4\tFALSE\tNA\tNA\tNA\tNA\tNA\tNA"};
  // Quality threshold 3-. Frequency threshold: 0.03. Total_DP = 3. Freq =
  // 0.666667 No Indel
  std::string t_20_03[3] = {
      "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0."
      "2\tFALSE\tgeneKatyPerry:id-test3\tAAG\tK\tATG\tM\t70",
      "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0.2\tFALSE\tid-"
      "testedit1\tGAA\tE\tGAT\tD\t70",
      "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0."
      "2\tFALSE\tgeneTaylorSwift:id-testedit2\tAGA\tR\tTGA\t*\t70"};
  // Quality threshold 25. Frequency threshold: 0.03. Total_DP = 2. Freq = 0.5
  // No Indel.
  std::string t_25_03[3] = {
      "test\t210\tA\tT\t1\t1\t41\t1\t1\t58\t0.5\t2\t0.4\tFALSE\tgeneKatyPerry:"
      "id-test3\tAAG\tK\tATG\tM\t70",
      "test\t210\tA\tT\t1\t1\t41\t1\t1\t58\t0.5\t2\t0.4\tFALSE\tid-"
      "testedit1\tGAA\tE\tGAT\tD\t70",
      "test\t210\tA\tT\t1\t1\t41\t1\t1\t58\t0.5\t2\t0."
      "4\tFALSE\tgeneTaylorSwift:id-testedit2\tAGA\tR\tTGA\t*\t70"};
  // Minimum depth threshold. Should be empty
  std::string t_25_03_20[] = {};
  num_success =
      call_var_check_outfile("../../data/test.indel", 20, 0, 0.02, t_20_02_1, 4);
  std::cout << num_success << std::endl;
  num_success +=
      call_var_check_outfile("../../data/test.indel2", 20, 0, 0.03, t_20_03, 3);
  std::cout << num_success << std::endl;
  num_success +=
      call_var_check_outfile("../../data/test.indel3", 25, 0, 0.03, t_25_03, 3);
  std::cout << num_success << std::endl;
  num_success += call_var_check_outfile("../../data/test.indel4", 25, 20, 0.03,
                                        t_25_03_20, 0);
  std::cout << num_success << std::endl;
  if (num_success == 0) return 0;
  return -1;
}
