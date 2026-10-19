#ifndef IVAR_SOLVE_CLUSTERING_H
#define IVAR_SOLVE_CLUSTERING_H
#include "gmm.h"

//Searches the effective GMM peak means for subsets that sum to 1, i.e. candidate
//population frequency solutions. Never touches variants.
class subset_sum_solver {
 private:
  std::vector<double> means;
  double unit_sum_error;
  std::vector<std::vector<double>> solution_sets;

  static bool within_error_range(const std::vector<double> &values, double target, double error);
  static double find_nearest_distance(const std::vector<double> &all_sums, double value);
  std::vector<std::vector<double>> find_solutions() const;
  std::vector<std::vector<double>> find_subsets(double target) const;
  bool account_for_inferred_means(const std::vector<double> &candidate) const;

 public:
  //the production unit sum tolerance, named so call sites don't hardcode it. Still
  //passed explicitly, since callers needing another value must be able to say so
  static constexpr double UNIT_SUM_ERROR = 0.10;
  //decomposition tolerance, the counterpart to UNIT_SUM_ERROR: how close a sum of
  //populations must land to a peak to explain it. Gates find_subsets,
  //account_for_inferred_means and variant_assigner::find_combination_peaks, so the
  //solver and the assigner always agree on which peaks a solution explains. Never
  //varied by a caller. TEST 8 exists to keep it distinct from UNIT_SUM_ERROR - do
  //not merge them.
  static constexpr double SUBSET_SUM_ERROR = 0.03;

  subset_sum_solver(std::vector<double> means, double unit_sum_error);
  bool solve();

  const std::vector<std::vector<double>> &get_solution_sets() const { return solution_sets; }
};

//Given a chosen solution, decides which consensus genomes each variant belongs to.
//Never does solution search.
class variant_assigner {
 private:
  std::vector<double> solution;
  std::vector<double> means;
  double threshold;

  std::vector<uint32_t> members_of(const std::vector<double> &subset) const;
  std::vector<std::vector<uint32_t>> find_combination_peaks(std::vector<double> &unresolved, std::vector<std::vector<uint32_t>> &ambiguous_indexes) const;
  void collect_matches(const std::vector<std::vector<uint32_t>> &inverse_groups, const std::vector<uint32_t> &top_clusters, std::vector<uint32_t> &out) const;

 public:
  variant_assigner(std::vector<double> solution, std::vector<double> means, double threshold);
  void assign(std::vector<variant> &variants) const;

  //remaps cluster_assigned out of model_means index space into eff_means space, a
  //preparation step for assign(). Static since it needs no solution.
  static void overwrite_cluster_assigned(std::vector<variant> &variants, const std::vector<double> &eff_means, const std::vector<double> &means);
};

#endif
