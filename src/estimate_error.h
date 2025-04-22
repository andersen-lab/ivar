#include <vector>
#include <fstream>
#include "./include/armadillo"
#include "gmm.h"
#ifndef estimate_error
#define estimate_error
double cluster_error(std::vector<variant> base_variants, uint8_t quality_threshold, uint32_t depth_cutoff);
#endif
