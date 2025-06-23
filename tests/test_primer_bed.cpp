#include <iostream>
#include <vector>

#include "../src/primer_bed.h"

int main() {
  int num_tests = 2;
  int success = 0;
  std::vector<primer>::iterator it;
  std::vector<primer> primers = populate_from_file("../data/test.bed");
  std::string primer_names[] = {"WNV_400_1_LEFT",     "WNV_400_1_LEFT_alt",
                                "WNV_400_2_LEFT",     "WNV_400_1_RIGHT",
                                "WNV_400_2_RIGHT",    "WNV_400_3_LEFT",
                                "WNV_400_2_LEFT_alt", "WNV_400_2_RIGHT_alt"};
  int primer_indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  unsigned int primer_start[] = {8, 7, 230, 359, 658, 569, 251, 352};
  unsigned int primer_end[] = {29, 27, 249, 380, 679, 590, 273, 378};
  int flag = 1;
  for (it = primers.begin(); it != primers.end(); ++it) {
    if (it->get_name() != primer_names[it - primers.begin()] ||
        it->get_start() != primer_start[it - primers.begin()] ||
        it->get_end() != primer_end[it - primers.begin()]) {
      std::cout << it->get_name() << std::endl;
      std::cout << it->get_end() << std::endl;
      flag = 0;
    }
  }
  for (it = primers.begin(); it != primers.end(); ++it) {
    if (it->get_indice() != primer_indices[it - primers.begin()]) {
      std::cout << "Wrong primer indice: " << it->get_name() << std::endl;
      flag = 0;
    }
  }
  success += flag;
  populate_pair_indices(primers, "../data/pair_information.tsv");
  int pair_indices[] = {3, -1, 4, 0, 2, -1, -1, -1};
  flag = 1;
  for (it = primers.begin(); it != primers.end(); ++it) {
    if (it->get_pair_indice() != pair_indices[it - primers.begin()]) {
      std::cout << "Wrong pair for " << it->get_name() << std::endl;
      flag = 0;
    }
  }
  success += flag;
  return (num_tests == success) ? 0 : -1;
}
