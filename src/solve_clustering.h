#include "gmm.h"
#ifndef solve_clustering
#define solve_clustering
void find_combinations(std::vector<double> means, uint32_t index, std::vector<double> &current, std::vector<std::vector<double>> &results, double error);
void solve_clusters(std::vector<variant> &variants, gaussian_mixture_model model, double estimated_error, std::vector<double> &solution);
std::vector<uint32_t> rewrite_amplicon_masking(std::vector<variant> variants, std::vector<double> means);
#endif
