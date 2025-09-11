#include <vector>
#include <fstream>
#include "./include/armadillo"
#include "gmm.h"
#ifndef estimate_error
#define estimate_error
std::vector<uint32_t>determine_outlier_points(std::vector<double> cluster, double threshold);
void cluster_error(std::vector<variant> base_variants, uint8_t quality_threshold, uint32_t depth_cutoff, double &error_rate);
#endif
