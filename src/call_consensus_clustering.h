#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

void call_majority_consensus(std::vector<variant> variants, uint32_t max_position, std::string clustering_file, double default_threshold);
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means, std::string ref);
#endif
