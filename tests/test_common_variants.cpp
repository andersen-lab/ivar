#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "../src/get_common_variants.h"

int main() {
  int num_success = 0;
  std::string files[] = {"../data/test.1.tsv", "../data/test.2.tsv",
                         "../data/test.3.tsv"};
  std::ifstream fin;
  std::map<std::string, std::string> file_tab_delimited_str;
  std::map<std::string, unsigned int> counts;
  // Read all files
  for (int i = 0; i < 3; ++i) {
    fin.open(files[i]);
    read_variant_file(fin, i, counts, file_tab_delimited_str);
    fin.close();
  }
  std::map<std::string, std::string>::iterator it =
      file_tab_delimited_str.begin();
  // Test tab delimited strings for test.1.tsv.
  std::string keys[] = {"test\t42\tG\tT\tid-test3\tAGG\tR\tATG\tM\t0",
                        "test\t69\tT\tG\tid-test3\tTTG\tL\tTGG\tW\t0",
                        "test\t320\tA\tT\tNA\tNA\tNA\tNA\tNA\t2"};
  std::string values[] = {"0\t0\t0\t1\t0\t49\t1\t1\t1\tFALSE",
                          "1\t0\t57\t1\t0\t53\t0.5\t2\t0.666667\tFALSE",
                          "1\t1\t35\t1\t1\t46\t0.5\t2\t0.666667\tFALSE"};
  int ctr = 0;
  while (it != file_tab_delimited_str.end() && ctr < 2) {
    // Check if key in map
    if (file_tab_delimited_str.find(keys[ctr]) ==
        file_tab_delimited_str.end()) {
      std::cout << keys[ctr] << std::endl;
      num_success = -1;
    } else {
      ctr++;
    }
    if ((it->first).compare(keys[ctr]) == 0 &&
        (it->second).compare(values[ctr]) !=
            0) {  // Check if value of key matches
      std::cout << it->first << ": " << it->second << std::endl;
      std::cout << keys[ctr] << ": " << values[ctr] << " -> Correct"
                << std::endl;
      num_success = -1;
    }
    it++;
  }
  // Test counts
  std::map<std::string, unsigned int>::iterator count_it = counts.begin();
  std::string count_keys[] = {"test\t42\tG\tT\tid-test3\tAGG\tR\tATG\tM\t",
                              "test\t320\tA\tT\tNA\tNA\tNA\tNA\tNA\t",
                              "test\t365\tA\tT\tNA\tNA\tNA\tNA\tNA\t",
                              "test\t42\tG\tT\tid-test4\tCAG\tQ\tCAT\tH\t"};
  unsigned int count_values[] = {1, 2, 3, 1};
  ctr = 0;
  while (count_it != counts.end() && ctr < 3) {
    if (counts.find(count_keys[ctr]) == counts.end()) {
      std::cout << count_keys[ctr] << std::endl;
      num_success = -1;
    } else {
      ctr++;
    }
    if ((count_it->first).compare(count_keys[ctr]) == 0 &&
        count_it->second != count_values[ctr]) {
      std::cout << count_it->first << ": " << count_it->second << std::endl;
      std::cout << count_keys[ctr] << ": " << count_values[ctr] << " -> Correct"
                << std::endl;
      num_success = -1;
    }
    count_it++;
  }
  if (num_success == 0) return 0;
  return -1;
}
