#include "gmm.h"
#ifndef solve_clustering
#define solve_clustering
double find_neighboring_cluster(double freq, uint32_t cluster_assigned, std::vector<double> means);
void find_combinations(std::vector<double> means, uint32_t index, std::vector<double> &current, std::vector<std::vector<double>> &results, double error);
bool test_cluster_deviation(double nearest_cluster, double variant_cluster, double std_dev);
void solve_clusters(std::vector<variant> &variants, gaussian_mixture_model model, double estimated_error, std::vector<double> &solution, std::string prefix, double default_threshold, uint32_t min_depth);
std::vector<uint32_t> rewrite_amplicon_masking(std::vector<variant> variants, std::vector<double> means);
void rewrite_position_masking(std::vector<variant> &variants);
void amplicon_specific_cluster_assignment(std::vector<variant> &variants, gaussian_mixture_model model);
void call_majority_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth);
std::string trim_trailing_ambiguities(std::string sequence, uint32_t max_position);
#endif
