#include <fstream>
#include <iostream>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "interval_tree.h"
#include "ref_seq.h"
#include "primer_bed.h"
#ifndef saga
#define saga

int preprocess_reads(std::string bam, std::string bed, std::string bam_out, std::string cmd, std::string pair_info, int32_t primer_offset, uint32_t min_depth, uint8_t min_qual);
double calculate_standard_deviation(std::vector<double> data);
#endif
