#include "gmm_1d.h"
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <numeric>
#include <random>
#include <iostream>
#include <algorithm>
#include <iterator>

constexpr double gmm_1d::PI;
constexpr double gmm_1d::DEFAULT_VAR_FLOOR;
constexpr double gmm_1d::MIN_BD_THRESHOLD;

double gmm_1d::log_normal_1d(double x, double mu, double var) {
  return -0.5 * (std::log(2 * PI * var) + std::pow(x - mu, 2) / var);
}

double gmm_1d::log_half_normal_1d(double x, double mu, double var, bool left_tail) {
  if(left_tail && x < mu)
    return -std::numeric_limits<double>::infinity();
  if (!left_tail && x > mu)
    return -std::numeric_limits<double>::infinity();

  return std::log(2.0) + log_normal_1d(x, mu, var);
}

// Using https://en.wikipedia.org/wiki/Digamma_function#Computation_and_approximation
double gmm_1d::digamma(double x) {
  double result = 0;
  while (x < 6) {
    result -= 1.0 / x;
    x += 1;
  }
  double inv_x  = 1.0 / x;
  double inv_x2 = inv_x * inv_x;
  result += std::log(x) - 0.5 * inv_x
            - inv_x2 * (1.0/12.0
            - inv_x2 * (1.0/120.0
            - inv_x2 * (1.0/252.0
            - inv_x2 * (1.0/240.0
            - inv_x2 * (1.0/132.0
            - inv_x2 * (691.0/32760.0))))));
  return result;
}

// Additional information at https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
double gmm_1d::log_sum_exp(const std::vector<double> &v) {
  double mx = *std::max_element(v.begin(), v.end());
  double s  = 0.0;
  for (double x : v) s += std::exp(x - mx);
  return mx + std::log(s);
}

void gmm_1d::initialize_k_means_1d(const std::vector<double> &x, int K, std::vector<int> &indices) {
  const size_t N = x.size();
  if (N == 0) {
    throw std::runtime_error("initialize_k_means_1d: no data available for k-means initialization");
  }
  if (K <= 0) {
    throw std::runtime_error("initialize_k_means_1d: number of indices must be positive");
  }
  if (static_cast<size_t>(K) > N) {
    throw std::runtime_error("initialize_k_means_1d: number of indices cannot exceed number of samples");
  }

  indices.clear();
  indices.reserve(K);

  std::uniform_int_distribution<int> uniform(0, N - 1);
  indices.push_back(uniform(rng));

  std::vector<double> dist(N);
  for (int c = 1; c < K; c++) {
    // Squared distance to nearest chosen center
    for (int i = 0; i < N; i++) {
      double min_d = std::numeric_limits<double>::max();
      for (int j : indices) {
        double d = (x[i] - x[j]) * (x[i] - x[j]);
        if (d < min_d) min_d = d;
      }
      dist[i] = min_d;
    }
    // Sample proportional to squared distance
    std::discrete_distribution<int> weighted(dist.begin(), dist.end());
    indices.push_back(weighted(rng));
  }
}

gmm_1d::Matrix gmm_1d::estimate_log_gaussian_prob(const std::vector<double>&x, const std::vector<double>& means, const std::vector<double>& inv_std_devs) const {
  int N = static_cast<int>(x.size());
  std::vector<double> prec(this->n_components), log_det(this->n_components);
  for (int k = 0; k < this->n_components; k++) {
    prec[k] = inv_std_devs[k] * inv_std_devs[k];
    log_det[k] = std::log(inv_std_devs[k]);
  }

  static const double LOG2PI = std::log(2.0 * gmm_1d::PI);
  Matrix result(N, std::vector<double>(this->n_components));
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < this->n_components; k++) {
      double diff = x[i] - means[k];
      double log_prob = prec[k] * diff * diff;
      result[i][k] = -0.5 * (LOG2PI + log_prob) + log_det[k];
    }
  }
  return result;
}

gmm_1d::Matrix gmm_1d::estimate_log_prob(const std::vector<double>&x) const {
  int N = static_cast<int>(x.size());
  Matrix log_gauss = estimate_log_gaussian_prob(x, means_, inv_std_devs_);

  std::vector<double> log_lambda(this->n_components);
  for (int k = 0; k < this->n_components; k++) {
    log_lambda[k] = std::log(2.0) + digamma(0.5 * degrees_of_freedom_[k]);
  }

  Matrix result(N, std::vector<double>(this->n_components));
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < this->n_components; k++) {
      result[i][k] = log_gauss[i][k] - 0.5 * std::log(degrees_of_freedom_[k]) + 0.5 * (log_lambda[k] - 1.0 / mean_precision_[k]);
    }
  }
  return result;
}

std::vector<double> gmm_1d::estimate_log_weights() const {
  std::vector<double> dg_sum(this->n_components), dg_a(this->n_components), dg_b(this->n_components);
  for (int k = 0; k < this->n_components; k++) {
    dg_sum[k] = digamma(this->stick_beta_a_[k] + this->stick_beta_b_[k]);
    dg_a[k] = digamma(this->stick_beta_a_[k]);
    dg_b[k] = digamma(this->stick_beta_b_[k]);
  }

  std::vector<double> prefix(this->n_components);
  double running = 0.0;
  for (int k = 0; k < this->n_components; k++) {
    prefix[k] = running;
    running += dg_b[k] - dg_sum[k];
  }

  std::vector<double> out(this->n_components);
  for (int k = 0; k < this->n_components; k++) {
    out[k] = dg_a[k] - dg_sum[k] + prefix[k];
  }
  return out;
}

gmm_1d::Matrix gmm_1d::estimate_weighted_log_prob(const std::vector<double>& X) const {
  Matrix lp = estimate_log_prob(X);
  std::vector<double> lw = estimate_log_weights();
  int N = static_cast<int>(X.size());
  Matrix out(N, std::vector<double>(this->n_components));
  for (int i = 0; i < N; i++)
    for (int k = 0; k < this->n_components; k++)
      out[i][k] = lp[i][k] + lw[k];
  return out;
}

void gmm_1d::e_step(const std::vector<double>& x, std::vector<double>& log_prob_norm, Matrix& log_resp) const {
  Matrix weighted = estimate_weighted_log_prob(x);
  int N = static_cast<int>(x.size());
  log_prob_norm.resize(N);
  log_resp.assign(N, std::vector<double>(this->n_components));
  for (int i = 0; i < N; i++) {
    log_prob_norm[i] = log_sum_exp(weighted[i]);
    for (int k = 0; k < this->n_components; k++)
      log_resp[i][k] = weighted[i][k] - log_prob_norm[i];
  }
}

void gmm_1d::estimate_gaussian_parameters(const std::vector<double>&x, const Matrix& resp, std::vector<double>& nk, std::vector<double>& means, std::vector<double>& covariances) const {
  int N = static_cast<int>(x.size());
  nk.assign(this->n_components, 0.0);
  means.assign(this->n_components, 0.0);
  covariances.assign(this->n_components, 0.0);
  std::vector<double> avg_X2(this->n_components, 0.0);

  for (int k = 0; k < this->n_components; k++) {
    for (int i = 0; i < N; i++) {
      nk[k] += resp[i][k];
      means[k] += resp[i][k] * x[i];
      avg_X2[k] += resp[i][k] * x[i] * x[i];
    }
    means[k] /= nk[k];
    avg_X2[k] /= nk[k];
    covariances[k] = avg_X2[k] - means[k] * means[k] + gmm_1d::REG_TERM;
  }
}

void gmm_1d::estimate_weights(const std::vector<double> &nk) {
  stick_beta_a_.resize(this->n_components);
  stick_beta_b_.resize(this->n_components);
  for (int k = 0; k < this->n_components; k++)
    stick_beta_a_[k] = 1.0 + nk[k];

  double suffix = 0.0;
  for (int k = this->n_components - 1; k >= 0; k--) {
    stick_beta_b_[k] = weight_concentration_prior_ + suffix;
    suffix += nk[k];
  }
}

void gmm_1d::estimate_means(const std::vector<double>& nk, const std::vector<double>& xk) {
  mean_precision_.resize(this->n_components);
  means_.resize(this->n_components);
  for (int k = 0; k < this->n_components; k++) {
    mean_precision_[k] = mean_precision_prior_ + nk[k];
    means_[k] = (mean_precision_prior_ * mean_prior_ + nk[k] * xk[k]) / mean_precision_[k];
  }
}

void gmm_1d::estimate_precisions(const std::vector<double>& nk,
                                 const std::vector<double>& xk,
                                 const std::vector<double>& sk) {
  degrees_of_freedom_.resize(this->n_components);
  variances_.resize(this->n_components);
  inv_std_devs_.resize(this->n_components);
  for (int k = 0; k < this->n_components; k++) {
    degrees_of_freedom_[k] = degrees_of_freedom_prior_ + nk[k];
    double diff = xk[k] - mean_prior_;
    variances_[k] = covariance_prior_ + nk[k] * (sk[k] + (mean_precision_prior_ / mean_precision_[k]) * diff * diff);
    variances_[k] /= degrees_of_freedom_[k];
    inv_std_devs_[k] = 1.0 / std::sqrt(variances_[k]);
  }
}

void gmm_1d::m_step(const std::vector<double>& x, const Matrix& log_resp) {
  int N = static_cast<int>(x.size());
  Matrix resp(N, std::vector<double>(this->n_components));
  for (int i = 0; i < N; i++)
    for (int k = 0; k < this->n_components; k++)
      resp[i][k] = std::exp(log_resp[i][k]);

  std::vector<double> nk, xk, sk;
  estimate_gaussian_parameters(x, resp, nk, xk, sk);
  estimate_weights(nk);
  estimate_means(nk, xk);
  estimate_precisions(nk, xk, sk);
}

double gmm_1d::compute_lower_bound(const Matrix& log_resp) const {
  int N = static_cast<int>(log_resp.size());

  // log_wishart_norm for n_features = 1
  double log_wishart = 0.0;
  for (int k = 0; k < this->n_components; k++) {
    double ld = std::log(inv_std_devs_[k]) - 0.5 * std::log(degrees_of_freedom_[k]);
    log_wishart += -(degrees_of_freedom_[k] * ld + degrees_of_freedom_[k] * 0.5 * std::log(2.0) + std::lgamma(0.5 * degrees_of_freedom_[k]));
  }

  double log_norm_weight = 0.0;
  for (int k = 0; k < this->n_components; k++)
    log_norm_weight -= std::lgamma(stick_beta_a_[k]) + std::lgamma(stick_beta_b_[k]) - std::lgamma(stick_beta_a_[k] + stick_beta_b_[k]);

  double entropy = 0.0;
  for (int i = 0; i < N; i++)
    for (int k = 0; k < this->n_components; k++)
      entropy -= std::exp(log_resp[i][k]) * log_resp[i][k];

  double sum_log_mp = 0.0;
  for (int k = 0; k < this->n_components; k++)
    sum_log_mp += std::log(mean_precision_[k]);

  return entropy - log_wishart - log_norm_weight - 0.5 * sum_log_mp;
}

void gmm_1d::initialize_parameters(const std::vector<double> &x) {
  if(x.empty()) {
    throw std::runtime_error("gmm_1d::initialize_parameters: input data is empty");
  }

  const int N = x.size();

  weight_concentration_prior_ = 1.0 / this->n_components;
  mean_precision_prior_ = 1.0;

  double s = 0.0;
  for (double x_ : x)
    s += x_;
  mean_prior_ = s / N;

  degrees_of_freedom_prior_ = 1.0;   // n_features = 1

  double var_sum = 0.0;
  for (double x_ : x) {
    double d = x_ - mean_prior_;
    var_sum += d * d;
  }
  covariance_prior_ = var_sum / (N - 1);
}

void gmm_1d::initialize(const std::vector<double>& x, const Matrix& resp) {
  std::vector<double> nk, means, covariances;
  estimate_gaussian_parameters(x, resp, nk, means, covariances);
  estimate_weights(nk);
  estimate_means(nk, means);
  estimate_precisions(nk, means, covariances);
}

void gmm_1d::compute_final_weights() {
  size_t K = stick_beta_a_.size();

  std::vector<double> weights(K);
  std::vector<double> tmp(K);
  std::vector<double> weight_dirichlet_sum(K);

  for (size_t k = 0; k < K; ++k) {
    weight_dirichlet_sum[k] = stick_beta_a_[k] + stick_beta_b_[k];
  }

  for (size_t k = 0; k < K; ++k) {
    tmp[k] = stick_beta_b_[k] / weight_dirichlet_sum[k];
  }

  std::vector<double> cumprod(K);
  cumprod[0] = 1.0;
  for (size_t k = 1; k < K; ++k) {
    cumprod[k] = cumprod[k - 1] * tmp[k - 1];
  }

  for (size_t k = 0; k < K; ++k) {
    weights[k] = (stick_beta_a_[k] / weight_dirichlet_sum[k]) * cumprod[k];
  }

  // normalize
  double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
  for (size_t k = 0; k < K; ++k) {
    weights[k] /= sum_weights;
  }

  this->weights_ = weights;
}

bool gmm_1d::fit(const std::vector<double>& x) {
  int N = x.size();
  initialize_parameters(x);

  std::vector<int> seed_indices;
  initialize_k_means_1d(x, this->n_components, seed_indices);

  Matrix resp(N, std::vector<double>(this->n_components, 0.0));
  for (int comp = 0; comp < this->n_components; comp++)
    resp[seed_indices[comp]][comp] = 1.0;

  initialize(x, resp);

  double prev_lower_bound = -std::numeric_limits<double>::infinity();
  std::vector<double> log_prob_norm;
  Matrix log_resp;

  for (int iter = 0; iter < gmm_1d::MAX_ITER; iter++) {
    e_step(x, log_prob_norm, log_resp);
    m_step(x, log_resp);
    double lb = compute_lower_bound(log_resp);
    this->elbo_history.push_back(lb);
    if (std::abs(lb - prev_lower_bound) <= gmm_1d::TOLERANCE) {
      this->converged = true;
      this->n_iter = iter + 1;
      break;
    }
    prev_lower_bound = lb;
  }

  if(!this->converged) {
    std::cerr << "[Warning]: Exceeded maximum iterations\n";
    this->n_iter = gmm_1d::MAX_ITER;
  }

  compute_final_weights();
  e_step(x, log_prob_norm, log_resp);

  labels_.clear();
  labels_.resize(N);
  for (int i = 0; i < N; i++) {
    int best = 0;
    for (int k = 1; k < this->n_components; k++)
      if (log_resp[i][k] > log_resp[i][best])
        best = k;
    labels_[i] = best;
  }

  return true;
}

std::vector<int> gmm_1d::predict(const std::vector<double>& x) const {
  Matrix w = estimate_weighted_log_prob(x);
  int N = static_cast<int>(x.size());
  std::vector<int> labels(N);
  for (int i = 0; i < N; i++) {
    int best = 0;
    for (int k = 1; k < this->n_components; k++)
      if (w[i][k] > w[i][best]) best = k;
    labels[i] = best;
  }
  return labels;
}

void gmm_1d::logit_transform(const std::vector<double> &x, std::vector<double> &transformed_x, double eps) {
  const size_t n = x.size();
  transformed_x.resize(n);

  for (size_t i = 0; i < n; ++i) {
    double v = x[i];
    v = std::min(std::max(v, eps), 1.0 - eps);
    transformed_x[i] = std::log(v / (1.0 - v));
  }
}

void gmm_1d::sigmoid_transform(const std::vector<double>& x, std::vector<double>& transformed_x, double eps) {
  transformed_x.resize(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    double s = 1.0 / (1.0 + std::exp(-x[i]));
    transformed_x[i] = std::min(std::max(s, eps), 1.0 - eps);
  }
}

// https://en.wikipedia.org/wiki/Bhattacharyya_distance
double gmm_1d::calculate_bhattacharyya_distance_1d(double mu1, double v1, double mu2, double v2) {
  v1 = std::max(v1, 1e-300);
  v2 = std::max(v2, 1e-300);
  const double term1 = 0.25 * std::log(0.25 * (v1 / v2 + v2 / v1 + 2.0));
  const double term2 = (mu1 - mu2) * (mu1 - mu2) / (4.0 * (v1 + v2));
  return term1 + term2;
}

// Based on Hennig et al. 2010
// https://doi.org/10.1007/s11634-010-0058-3
int gmm_1d::get_distinct_components_count(const std::vector<uint32_t>& sites, double min_bd_threshold) {
  const int G = n_components;
  if (G <= 0) return 0;

  // Get minimum number of components from site constraint
  int m_max = 0;
  {
    std::vector<uint32_t> s = sites;
    std::sort(s.begin(), s.end());
    int run = 0;
    for (size_t i = 0; i < s.size(); ++i) {
      if (i == 0 || s[i] == s[i - 1]) {
        run++;
      } else {
        m_max = std::max(m_max, run);
        run = 1;
      }
    }
    if (!s.empty())
      m_max = std::max(m_max, run);
  }

  std::vector<std::vector<int>> clusters(G);
  for (int g = 0; g < G; ++g)
    clusters[g].push_back(g);

  // BD distance matrix
  std::vector<std::vector<double>> dist_matrix(G, std::vector<double>(G));
  for (int i = 0; i < G; ++i) {
    for (int j = 0; j < G; ++j) {
      if (i == j) {
        dist_matrix[i][j] = 0.0;
      } else {
        dist_matrix[i][j] = calculate_bhattacharyya_distance_1d(
            means_[i], variances_[i], means_[j], variances_[j]
        );
      }
    }
  }

  // Agglomerative clustering but merge closest clusters first
  while (clusters.size() > 1) {
    double min_dist = std::numeric_limits<double>::infinity();
    int merge_i = -1, merge_j = -1;

    for (int i = 0; i < clusters.size(); ++i) {
      for (int j = i + 1; j < clusters.size(); ++j) {
        double cluster_dist = std::numeric_limits<double>::infinity();
        for (int ci : clusters[i]) {
          for (int cj : clusters[j]) {
            cluster_dist = std::min(cluster_dist, dist_matrix[ci][cj]);
          }
        }

        if (cluster_dist < min_dist) {
          min_dist = cluster_dist;
          merge_i = i;
          merge_j = j;
        }
      }
    }

    // Stop if minimum distance across all clusters exceeds threshold
    if (min_dist >= min_bd_threshold) {
      break;
    }

    // Merge cluster j into cluster i
    clusters[merge_i].insert(
        clusters[merge_i].end(),
        clusters[merge_j].begin(),
        clusters[merge_j].end()
    );

    clusters.erase(clusters.begin() + merge_j);
  }

  std::vector<double> m_sigmoid(means_.size());
  sigmoid_transform(means_, m_sigmoid);

  for (int c = 0; c < clusters.size(); ++c) {

    int rep = clusters[c][0];
    double max_w = weights_[rep];

    for (int g : clusters[c]) {
      if (weights_[g] > max_w) {
        max_w = weights_[g];
        rep = g;
      }
    }
  }

  return clusters.size();
}

std::vector<int> gmm_1d::get_effective_components(std::vector<int> labels) const {
  std::vector<int> tmp = labels;
  std::sort(tmp.begin(), tmp.end());

  std::vector<int> result;
  result.reserve(tmp.size());

  for (double x : tmp) {
    if (result.empty() || x != result.back()) {
      result.push_back(x);
    }
  }
  return result;
}

std::vector<double> gmm_1d::get_effective_means(std::vector<int> component_indices) const {
  std::vector<double> eff_means;
  for (int i : component_indices) {
    eff_means.push_back(this->means_[i]);
  }
  return eff_means;
}

std::vector<double> gmm_1d::get_effective_vars(std::vector<int> component_indices) const {
  std::vector<double> eff_variances;
  for (int i : component_indices) {
    eff_variances.push_back(this->variances_[i]);
  }
  return eff_variances;
}

std::vector<double> gmm_1d::get_effective_weights(std::vector<int> component_indices) const {
  std::vector<double> eff_weights;
  for (int i : component_indices) {
    eff_weights.push_back(this->weights_[i]);
  }
  return eff_weights;
}
