#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

void cluster_consensus(std::vector<variant> variants, std::string clustering_file, std::string variants_filename, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means);
#endif
