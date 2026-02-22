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
  static constexpr double MIN_BD_THRESHOLD = 1.0; // -log(d*) roughly based on Hennig et al. 2010s
  static constexpr double TOLERANCE = 1e-6;
  static constexpr double REG_TERM = 1e-6;
  static constexpr int MAX_ITER = 200;

  int n_components;
  unsigned int seed;
  std::mt19937 rng;

  // Priors
  double weight_concentration_prior_;
  double mean_precision_prior_;
  double mean_prior_;
  double degrees_of_freedom_prior_;
  double covariance_prior_;

  // Model parameters
  std::vector<double> stick_beta_a_, stick_beta_b_;
  std::vector<double> mean_precision_;
  std::vector<double> means_;
  std::vector<double> degrees_of_freedom_;
  std::vector<double> variances_;
  std::vector<double> inv_std_devs_;

  // Computed parameters
  std::vector<double> weights_;
  std::vector<double> elbo_history;
  bool converged = false;
  int n_iter = 0;
  std::vector<int> labels_;

  // Helpers
  static double digamma(double x);
  static double log_sum_exp(const std::vector<double> &v);
  static double log_normal_1d(double x, double mu, double var);
  static double log_half_normal_1d(double x, double mu, double var, bool left_tail);
  void initialize_k_means_1d(const std::vector<double>& x, int K, std::vector<int>& indices);

  // E-step
  using Matrix = std::vector<std::vector<double>>;
  Matrix estimate_log_gaussian_prob(const std::vector<double>& x, const std::vector<double>& means, const std::vector<double>& inv_std_devs) const;
  Matrix estimate_log_prob(const std::vector<double>& x) const;
  std::vector<double> estimate_log_weights() const;
  Matrix estimate_weighted_log_prob(const std::vector<double>& x) const;
  void e_step(const std::vector<double>& x, std::vector<double>& log_prob_norm, Matrix& log_resp) const;

  // M-step components
  void estimate_gaussian_parameters(const std::vector<double>& x, const Matrix& resp, std::vector<double>& nk, std::vector<double>& means, std::vector<double>& variances) const;
  void estimate_weights(const std::vector<double>& nk);
  void estimate_means(const std::vector<double>& nk, const std::vector<double>& xk);
  void estimate_precisions(const std::vector<double>& nk, const std::vector<double>& xk, const std::vector<double>& sk);
  void m_step(const std::vector<double>& x, const Matrix& log_resp);

  double compute_lower_bound(const Matrix& log_resp) const;
  void initialize_parameters(const std::vector<double>& x);
  void initialize(const std::vector<double>& x, const Matrix& resp);
  void compute_final_weights();

 public:

  explicit gmm_1d(int n_components = 2, unsigned int seed = std::random_device{}()) : seed(seed), rng(seed), n_components(n_components) {}
  bool fit(const std::vector<double>& x);
  std::vector<int> predict(const std::vector<double> &x) const;

  static void logit_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps = 1e-6);
  static void sigmoid_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps = 1e-6);
  static double calculate_bhattacharyya_distance_1d(double mu1, double v1, double mu2, double v2);
  int get_distinct_components_count(const std::vector<uint32_t>& sites, double min_bd_threshold);

  std::vector<double> get_weights() const { return weights_; }
  std::vector<double> get_means() const { return means_; }
  std::vector<double> get_variances() const { return variances_; }

  std::vector<int> get_effective_components(std::vector<int> labels) const;
  std::vector<double> get_effective_means(std::vector<int> labels) const;
  std::vector<double> get_effective_vars(std::vector<int> labels) const;
  std::vector<double> get_effective_weights(std::vector<int> labels) const;

  std::vector<double> get_elbo_history() const { return elbo_history; }
  bool is_converged() const { return converged; }
  int get_n_iter() const { return n_iter; }

  void set_seed(unsigned int s) {
      seed = s;
      rng.seed(seed);
  }

  unsigned int get_seed() const {
      return seed;
  }
};

#endif
