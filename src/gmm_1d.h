#ifndef IVAR_GMM_1D_H
#define IVAR_GMM_1D_H
#include <cstdint>
#include <cmath>
#include <vector>

// As described in Sanderson et al. 2017
// https://arma.sourceforge.net/armadillo_spcs_2017.pdf

class gmm_1d {
 private:
  // From https://learn.microsoft.com/en-us/cpp/c-runtime-library/math-constants
  static constexpr double PI =	3.14159265358979323846;
  static constexpr double DEFAULT_VAR_FLOOR = 1e-3;
  static constexpr double DEFAULT_WEIGHT_FLOOR = 1e-3;
  static constexpr double MIN_BD_THRESHOLD = 5.3025850929940455; // -log(d*) roughly based on Hennig et al. 2010s

  static double log_normal_1d(double x, double mu, double var);

  static double log_sum_exp(const double* x, size_t n);
  static double log_sum_exp(std::vector<double> x) {
    return log_sum_exp(x.data(), x.size());
  }


  static void generate_assignments(int G, int m, std::vector<std::vector<int>>& assignments);

  std::vector<double> means;
  std::vector<double> vars;
  std::vector<double> weights;
  int n_components;
  double var_floor = DEFAULT_VAR_FLOOR;
  double weight_floor = DEFAULT_WEIGHT_FLOOR;

  std::pair<std::vector<std::vector<double>>, double> site_resp_constrained_by_site(const std::vector<std::vector<double>>& logA) const;

  double e_step_1d(const std::vector<double>& x, const std::vector<uint32_t>& site_id, std::vector<std::vector<double>>& resp) const;
  void m_step_1d(const std::vector<double>& x, const std::vector<std::vector<double>>& resp);
  void initialize_k_means_1d(const std::vector<double>& x, uint32_t K, bool adaptive, uint32_t n_iter = 10, uint32_t seed = 112358);

 public:

  explicit gmm_1d(uint32_t n_components = 2) : n_components(n_components) {}
  static void logit_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps = 1e-6);
  static void sigmoid_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps = 1e-6);
  static double calculate_bhattacharyya_distance_1d(double mu1, double v1, double mu2, double v2);
  uint32_t get_distinct_components_count(const std::vector<uint32_t>& sites, double min_bd_threshold = MIN_BD_THRESHOLD);

  int fit(const std::vector<double>& x, const std::vector<uint32_t>& sites, std::vector<double>& logL_history, int n_iter = 20, double tolerance = -1.0, bool adaptive = false, bool logging = true, unsigned int seed = 112358);
  bool predict(const std::vector<double>& x, const std::vector<uint32_t>& sites, std::vector<uint32_t>& assigned_components, std::vector<std::vector<double>>& marginal_posterior_probabilities) const;

  std::vector<double> get_weights() const { return weights; }
  std::vector<double> get_means() const { return means; }
  std::vector<double> get_vars() const { return vars; }
  double get_log_likelihood(const std::vector<double>& x, const std::vector<uint32_t>& sites) const;
  double get_bic(const std::vector<double>& x, const std::vector<uint32_t>& sites) const;

  void set_var_floor(double val) { var_floor = val; }
  void set_weight_floor(double val) { weight_floor = val; }

  //store the cluster infor post-merging
  std::vector<double> merged_means;
  std::vector<double> merged_vars;
  std::vector<double> merged_weights; 
};

#endif
