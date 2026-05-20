#include "gmm.h"
#ifndef solve_clustering
#define solve_clustering

void overwrite_cluster_assigned(std::vector<variant> &variants, 
                                std::vector<double> eff_means, 
                                std::vector<double> means);
bool subset_sum(std::vector<double> &means, std::vector<std::vector<double>> &solution_sets, const double error);
double find_neighboring_cluster(double freq, uint32_t cluster_assigned, std::vector<double> means);
void find_combinations(std::vector<double> means, uint32_t index, std::vector<double> &current, std::vector<std::vector<double>> &results, double error);
bool test_cluster_deviation(double nearest_cluster, double variant_cluster, double std_dev);
std::vector<uint32_t> rewrite_amplicon_masking(std::vector<variant> variants, std::vector<double> means);
void rewrite_position_masking(std::vector<variant> &variants);
void call_majority_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth);
void assign_variants_solution(std::vector<double> solution, std::vector<variant> &variants, std::vector<double> means);
std::vector<std::vector<double>> subset_sum(std::vector<double> means);
#endif