#include <vector>
#include <fstream>
#include "./include/armadillo"
#ifndef estimate_error
#define estimate_error
double cluster_error(std::string filename, double quality_threshold, uint32_t depth_cutoff);
#endif
