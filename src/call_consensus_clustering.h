#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

void rewrite_position_masking(std::vector<variant> variants);
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, std::string variants_filename, double default_threshold, uint32_t min_depth, uint8_t min_qual);
void find_combinations(std::vector<float> means, uint32_t index, std::vector<float> &current, std::vector<std::vector<float>> &results, float error);
#endif
