#include<iostream>
#include "../src/call_consensus_pileup.h"
#include "../src/allele_functions.h"

int main() {
  int num_tests = 10;
  allele a1 = {
    "T",
    8,
    0,
    30,
    0,
    0,
    0
  };
  allele ai1 = {
    "+A",
    6,
    0,
    30,
    0,
    0,
    0
  };
  allele a2 = {
    "A",
    2,
    0,
    30,
    0,
    0,
    0
  };
  int success = 0;
  allele arr[] = {a1, ai1, a2};
  ret_t s;
  std::vector<allele> ad(arr, arr+sizeof(arr)/sizeof(allele));
  int size = sizeof(arr)/sizeof(allele);
  for(int i = 0;i<size;i++){
    ad.at(i) = arr[i];
  }

  s = get_consensus_allele(ad,20,.8, 'N', 0, 0);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("TN") == 0) ? 1: 0;
  success += (s.q.compare("??") == 0) ? 1 : 0;

  s = get_consensus_allele(ad,20,.8, 'N', 1, 1);
  std::cout << s.nuc << ": " << s.q << std::endl; 
  success += (s.nuc.compare("T") == 0) ? 1: 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;
 
  s = get_consensus_allele(ad,20,.8, 'N',.6, 0.6);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("TN") == 0) ? 1: 0;
  success += (s.q.compare("??") == 0) ? 1 : 0;
  
  s = get_consensus_allele(ad,20,.8, 'N',.8, 0.8);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("T") == 0) ? 1: 0;
  success += (s.q.compare("?") == 0) ? 1 : 0;
    
  //adding a test to check snv against insertion threshold
  allele a3 = {
    "T",
    5,
    0,
    30,
    0,
    0,
    0
  };
  allele a4 = {
    "+GGGGGGG",
    3,
    0,
    30,
    0,
    0,
    0
  };
  allele a5 = {
    "A",
    2,
    0,
    30,
    0,
    0,
    0
  };

  allele arr1[] = {a3, a4, a5};
  std::vector<allele> ad1(arr1, arr1+sizeof(arr1)/sizeof(allele));
  
 for(int i = 0;i<size;i++){
    ad1.at(i) = arr1[i];
  }

  //not bothering with single insertion, should be checked elsewhere
  //first param is 1bp second is >1bp 
  s = get_consensus_allele(ad1,20,0, 'N',0.8, 0);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("TGGGGGGG") == 0) ? 1: 0;
   
  s = get_consensus_allele(ad1,20, 0, 'N',0, 0.8);
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("T") == 0) ? 1: 0;
  
  return (success == num_tests) ? 0 : -1;
}
