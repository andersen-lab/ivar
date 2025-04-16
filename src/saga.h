#include <fstream>
#include <iostream>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "interval_tree.h"
#include "primer_bed.h"
#include "trim_primer_quality.h"
#ifndef saga
#define saga

double calculate_standard_deviation(std::vector<double> data);
int preprocess_reads(std::string bam, std::string bed, std::string bam_out, std::string cmd, std::string pair_info, int32_t primer_offset, uint32_t min_depth, uint8_t min_qual);
double calculate_standard_deviation(std::vector<double> data);
#endif
