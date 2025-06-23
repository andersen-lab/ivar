#include <libgen.h>
#include <stdint.h>
#include <string.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "allele_functions.h"

#ifndef call_consensus_from_pileup
#define call_consensus_from_pileup

struct ret_t {
  std::string nuc;
  std::string q;
};

void format_alleles(std::vector<allele> &ad);
int call_consensus_from_plup(std::istream &cin, std::string seq_id,
                             std::string out_file, uint8_t min_qual,
                             double threshold, uint32_t min_depth, char gap,
                             bool min_coverage_flag,
                             double min_insert_threshold);
ret_t get_consensus_allele(std::vector<allele> ad, uint8_t min_qual,
                           double threshold, char gap,
                           double min_insert_threshold);

#endif
