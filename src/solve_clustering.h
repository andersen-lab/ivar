#include "gmm.h"
#ifndef solve_clustering
#define solve_clustering

void overwrite_cluster_assigned(std::vector<variant> &variants,
                                std::vector<double> eff_means,
                                std::vector<double> means);
bool subset_sum(std::vector<double> &means, std::vector<std::vector<double>> &solution_sets, const double error);
void find_combinations(std::vector<double> means, uint32_t index, std::vector<double> &current, std::vector<std::vector<double>> &results, double error);
void rewrite_position_masking(std::vector<variant> &variants);
void assign_variants_solution(std::vector<double> solution, std::vector<variant> &variants, std::vector<double> means, double threshold);
#endif