#include <iostream>

#include "../src/allele_functions.h"
#include "../src/call_consensus_pileup.h"

int main() {
  int num_tests = 8;
  allele a1 = {"A", 40, 0, 30, 0, 0, 0};
  allele a2 = {"T", 40, 10, 30, 0, 0, 0};
  allele a3 = {"G", 40, 10, 10, 0, 0, 0};
  allele a4 = {"+T", 30, 10, 10, 0, 0, 0};
  allele a5 = {"+AG", 30, 10, 10, 0, 0, 0};
  allele a6 = {"+AC", 30, 4, 10, 0, 0, 0};
  allele a7 = {"+WT", 30, 2, 10, 0, 0, 0};
  int success = 0;
  allele arr[] = {a1, a2, a3, a4, a5};
  ret_t s;
  std::vector<allele> ad(arr, arr + sizeof(arr) / sizeof(allele));
  print_allele_depths(ad);
  int size = sizeof(arr) / sizeof(allele);
  for (int i = 0; i < size; i++) {
    ad.at(i) = arr[i];
  }
  s = get_consensus_allele(ad, 0, 0, 'N', 0.5);
  success += (s.nuc.compare("DW") == 0) ? 1 : 0;
  success += (s.q.compare("8+") == 0) ? 1 : 0;
  ad.push_back(a6);
  s = get_consensus_allele(ad, 0, 0, 'N', 0);
  success += (s.nuc.compare("DAS") == 0) ? 1 : 0;
  success += (s.q.compare("8++") == 0) ? 1 : 0;
  ad.push_back(a7);
  s = get_consensus_allele(ad, 0, 0, 'N', 0);
  success += (s.nuc.compare("DAB") == 0) ? 1 : 0;
  success += (s.q.compare("8++") == 0) ? 1 : 0;
  allele a8 = {"-AT", 10, 5, 20, 0, 0, 0};
  allele a9 = {"-A", 10, 5, 20, 0, 0, 0};
  allele a10 = {"-AG", 10, 5, 20, 0, 0, 0};
  allele a11 = {"-AGT", 10, 5, 20, 0, 0, 0};

  allele del_arr[] = {a8, a9, a10, a11};
  ret_t del_s;
  std::vector<allele> del_ad(del_arr,
                             del_arr + sizeof(del_arr) / sizeof(allele));
  s = get_consensus_allele(del_ad, 0, 0, 'N', 1);
  success += (s.nuc.compare("") == 0) ? 1 : 0;
  success += (s.q.compare("") == 0) ? 1 : 0;
  return (success == num_tests) ? 0 : -1;
}
