#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

std::string trim_leading_ambiguities(std::string sequence, uint32_t min_position);
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means, std::string ref, double error_rate);
#endif
