#include "allele_functions.h"
#include "interval_tree.h"
#ifndef genomic_position
#define genomic_position

class genomic_position {
  public:
    struct amplicon_info {
      ITNode* node;
      std::vector<allele> amp_alleles;
      uint32_t amp_depth;
    };
    std::vector<amplicon_info> amplicons;
    bool overlap = false;
    uint32_t pos;
    uint32_t depth=0;
    std::vector<allele> alleles;
    std::unordered_map<std::string, std::vector<double>> amplicon_frequencies = {};
    std::vector<uint32_t> amplicon_numbers; //amplicons that cover this position
    void update_alleles(std::string nt, uint32_t count, uint32_t qual);
    bool flux;
    bool in_primer=false;
};
int check_position_exists(uint32_t p, std::vector<position> positions);

#endif
