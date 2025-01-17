#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

void cluster_consensus(std::vector<variant> variants, std::string clustering_file, std::string variants_filename);
void find_combinations(std::vector<float> means, uint32_t index, std::vector<float> &current, std::vector<std::vector<float>> &results);
#endif
