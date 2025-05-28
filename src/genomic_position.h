#ifndef genomic_position_h
#define genomic_position_h
#include <vector>
#include <string>
#include "allele_functions.h"

class ITNode; //forward delcaration
struct amplicon_info {
  ITNode* node = nullptr;
  std::vector<allele> amp_alleles;
  uint32_t amp_depth = 0;
  void update_alleles(std::string allele, uint32_t qual);
};

class genomic_position {
  public:
    std::vector<amplicon_info> amplicons;
    uint32_t pos;
    uint32_t depth=0;
    std::vector<allele> alleles;
    std::unordered_map<std::string, std::vector<double>> amplicon_frequencies = {};
    std::vector<uint32_t> amplicon_numbers; //amplicons that cover this position
    void update_alleles(std::string nt, uint32_t qual);
    bool flux=false; //position frequency may fluctuate between amplicons
    bool overlap=false; //positon is covered by multiple amplicons

};
void combine_haplotypes(std::vector<genomic_position> &global_positions);
void populate_positions(std::vector<genomic_position> &positions, uint32_t max_position);
int check_position_exists(uint32_t p, std::vector<genomic_position> positions);
void assign_read(ITNode *node, std::vector<uint32_t> final_positions, std::vector<std::string> final_bases, std::vector<uint32_t> final_qualities, std::vector<genomic_position> &global_positions);
void add_variants(std::vector<uint32_t> final_positions, std::vector<std::string> final_bases, std::vector<uint32_t> final_qualities, std::vector<genomic_position> &global_positions);
#endif
