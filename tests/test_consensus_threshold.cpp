#include <iostream>

#include "../src/allele_functions.h"
#include "../src/call_consensus_pileup.h"

int main() {
  int num_tests = 6;
  allele a1 = {"A", 6, 0, 30, 0, 0, 0};
  allele ai1 = {"+T", 3, 0, 30, 0, 0, 0};
  allele a2 = {"T", 1, 10, 30, 0, 0, 0};
  allele a3 = {"C", 1, 10, 30, 0, 0, 0};
  allele a4 = {"G", 1, 10, 30, 0, 0, 0};
  allele a5 = {"*", 1, 10, 30, 0, 0, 0};
  allele a6 = {"G", 2, 10, 30, 0, 0, 0};
  int success = 0;
  allele arr[] = {a1, ai1, a2, a3, a4, a5};
  ret_t s;
  std::vector<allele> ad(arr, arr + sizeof(arr) / sizeof(allele));
  int size = sizeof(arr) / sizeof(allele);
  for (int i = 0; i < size; i++) {
    ad.at(i) = arr[i];
  }
  s = get_consensus_allele(ad, 20, .6, 'N', 0.6);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("A") == 0) ? 1 : 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;
  s = get_consensus_allele(ad, 20, .7, 'N', 1);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("N") == 0) ? 1 : 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;
  ad.erase(ad.begin() + 4, ad.begin() + 5);
  ad.push_back(a6);
  s = get_consensus_allele(ad, 20, .7, 'N', 1);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("R") == 0) ? 1 : 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;
  return (success == num_tests) ? 0 : -1;
}
