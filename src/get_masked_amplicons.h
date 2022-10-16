#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "primer_bed.h"

#ifndef get_masked_amplicons
#define get_masked_amplicons

int get_primers_with_mismatches(std::string bed, std::string vpath, std::string out, std::string primer_pair_file);
#endif
