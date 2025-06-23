#include <iostream>

#include "../src/allele_functions.h"
#include "../src/call_consensus_pileup.h"

int main() {
  int num_tests = 8;
  allele a1 = {"T", 8, 0, 30, 0, 0, 0};
  allele ai1 = {"+A", 6, 0, 30, 0, 0, 0};
  allele a2 = {"A", 2, 0, 30, 0, 0, 0};
  int success = 0;
  allele arr[] = {a1, ai1, a2};
  ret_t s;
  std::vector<allele> ad(arr, arr + sizeof(arr) / sizeof(allele));
  int size = sizeof(arr) / sizeof(allele);
  for (int i = 0; i < size; i++) {
    ad.at(i) = arr[i];
  }

  s = get_consensus_allele(ad, 20, .8, 'N', 0);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("TN") == 0) ? 1 : 0;
  success += (s.q.compare("??") == 0) ? 1 : 0;

  s = get_consensus_allele(ad, 20, .8, 'N', 1);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("T") == 0) ? 1 : 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;

  s = get_consensus_allele(ad, 20, .8, 'N', .6);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("TN") == 0) ? 1 : 0;
  success += (s.q.compare("??") == 0) ? 1 : 0;

  s = get_consensus_allele(ad, 20, .8, 'N', .8);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("T") == 0) ? 1 : 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;

  return (success == num_tests) ? 0 : -1;
}
