#include <vector>
#include <fstream>
#include "./include/armadillo"
#ifndef estimate_error
#define estimate_error
std::vector<float> cluster_error(std::string filename, double quality_threshold);
#endif
