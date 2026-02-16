#ifndef IVAR_GMM_1D_H
#define IVAR_GMM_1D_H
#include <cmath>
#include <vector>
#include <random>

// Based on Sanderson et al. 2017
// https://arma.sourceforge.net/armadillo_spcs_2017.pdf

class gmm_1d {
 private:
  // From https://learn.microsoft.com/en-us/cpp/c-runtime-library/math-constants
  static constexpr double PI =	3.14159265358979323846;
  static constexpr double DEFAULT_VAR_FLOOR = 1e-3;
  static constexpr double DEFAULT_WEIGHT_FLOOR = 1e-3;
  static constexpr double MIN_BD_THRESHOLD = 9.3; // -log(d*) roughly based on Hennig et al. 2010s

  static double log_normal_1d(double x, double mu, double var);
  static double log_half_normal_1d(double x, double mu, double var, bool left_tail);

  static double log_sum_exp(const double* x, size_t n);
  static double log_sum_exp(std::vector<double> x) {
    return log_sum_exp(x.data(), x.size());
  }


  static void generate_assignments(int G, int m, std::vector<std::vector<int>>& assignments);
  static void backtrack_assignments(int G, int m, std::vector<std::vector<int>>& assignments, std::vector<int>& current, std::vector<bool>& used);

  std::vector<double> means;
  std::vector<double> vars;
  std::vector<double> weights;
  std::vector<double> data_weights;
  unsigned int seed;
  std::mt19937 rng;

  int n_components;
  double var_floor = DEFAULT_VAR_FLOOR;
  double weight_floor = DEFAULT_WEIGHT_FLOOR;
  bool use_half_normal_for_noise = true;

  double HALF_NORMAL_LEFT_THRESHOLD = 0.01;
  double HALF_NORMAL_RIGHT_THRESHOLD = 0.99;


  enum ComponentType { GAUSSIAN, HALF_NORMAL_LEFT, HALF_NORMAL_RIGHT };
  std::vector<ComponentType> component_types;

  std::pair<std::vector<std::vector<double>>, double> site_resp_constrained_by_site(const std::vector<std::vector<double>>& logA, const std::vector<double>& site_weights) const;

  double e_step_1d(const std::vector<double>& x, const std::vector<uint32_t>& site_id, std::vector<std::vector<double>>& resp) const;
  void m_step_1d(const std::vector<double>& x, const std::vector<std::vector<double>>& resp);
  void initialize_k_means_1d(const std::vector<double>& x_filtered, int K, std::vector<double>& centers, int n_iter = 10);
  void compute_data_weights(const std::vector<uint32_t>& depths);

 public:

  explicit gmm_1d(int n_components = 2, unsigned int seed = std::random_device{}()) : seed(seed), rng(seed), n_components(n_components) {}
  static void logit_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps = 1e-6);
  static void sigmoid_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps = 1e-6);
  static double calculate_bhattacharyya_distance_1d(double mu1, double v1, double mu2, double v2);
  int get_distinct_components_count(const std::vector<uint32_t>& sites, double min_bd_threshold = MIN_BD_THRESHOLD);

  bool fit(const std::vector<double>& x, const std::vector<uint32_t>& sites, std::vector<double>& logL_history, const std::vector<uint32_t>& depths = {}, int n_iter = 20, double tolerance = -1.0, bool adaptive = false);
  bool predict(const std::vector<double>& x, const std::vector<uint32_t>& sites, std::vector<int>& assigned_components, std::vector<std::vector<double>>& marginal_posterior_probabilities) const;

  std::vector<double> get_weights() const { return weights; }
  std::vector<double> get_means() const { return means; }
  std::vector<double> get_vars() const { return vars; }
  double get_log_likelihood(const std::vector<double>& x, const std::vector<uint32_t>& sites) const;
  double get_bic(const std::vector<double>& x, const std::vector<uint32_t>& sites) const;

  void set_var_floor(double val) { var_floor = val; }
  void set_weight_floor(double val) { weight_floor = val; }
  void set_use_half_normal_for_noise(bool val) { use_half_normal_for_noise = val; }
  void set_half_normal_thresholds(double left_threshold, double right_threshold) {
      HALF_NORMAL_LEFT_THRESHOLD = left_threshold;
      HALF_NORMAL_RIGHT_THRESHOLD = right_threshold;
  }
  void set_seed(unsigned int s) {
      seed = s;
      rng.seed(seed);
  }
  unsigned int get_seed() const {
      return seed;
  }

  //store the cluster infor post-merging
  std::vector<double> merged_means;
  std::vector<double> merged_vars;
  std::vector<double> merged_weights;
};

#endif
