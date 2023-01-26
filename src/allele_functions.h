#include<string>
#include<vector>
#include<stdint.h>
#include<unordered_map>
#ifndef allele_functions
#define allele_functions

struct allele{
  std::string nuc;
  double depth;
  uint32_t reverse;
  uint8_t mean_qual;
  uint32_t beg;
  uint32_t end;
  float tmp_mean_qual;
  bool ref = false; //is this the reference  
  bool operator < (const allele& a) const{
    return (nuc.compare(a.nuc) > 0) ? true : false;
  }
  bool operator == (const allele& a) const{
    return (nuc.compare(a.nuc) == 0) ? true : false;
  }
};

//data struct to hold all alleles at a position
struct position{
  std::vector<allele> ad; //all possible alleles at this position
  uint32_t pos = 0; //positions relative to the reference
  double depth = 0; //total ungapped depth at the position
};

int check_allele_exists(std::string n, std::vector<allele> ad);
std::vector<allele> update_allele_depth(char ref,std::string bases, std::string qualities, uint8_t min_qual);
//for use in consensus thresholding
void update_allele_depth(std::vector<position> &all_positions, std::vector<int> haplotypes, std::vector<uint32_t> positions, std::vector<float> qualities, std::unordered_map<int, std::string> dict_decode);
void print_allele_depths(std::vector<allele> ad);
void print_single_allele(allele a);
int find_ref_in_allele(std::vector<allele> ad, char ref);
char gt2iupac(char a, char b);
char codon2aa(char n1, char n2, char n3);

#endif
