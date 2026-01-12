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
constexpr double gmm_1d::DEFAULT_WEIGHT_FLOOR;
constexpr double gmm_1d::MIN_BD_THRESHOLD;

double gmm_1d::log_normal_1d(double x, double mu, double var) {
  return -0.5 * (std::log(2 * PI * var) + std::pow(x - mu, 2) / var);
}

// Additional information at https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
double gmm_1d::log_sum_exp(const double *x, size_t n) {
  double max_val = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    if (x[i] > max_val) {
      max_val = x[i];
    }
  }

  if (!std::isfinite(max_val))
    return max_val;

  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += std::exp(x[i] - max_val);
  }

  return max_val + std::log(sum);
}

// Assign G components to m variants per site under the constraint that m < G
void gmm_1d::generate_assignments(int G, int m, std::vector<std::vector<int>> &assignments) {
  std::vector<int> perm(G);
  for (int i = 0; i < G; ++i) perm[i] = i;
  do {
    assignments.emplace_back(perm.begin(), perm.begin() + m);
  } while (std::next_permutation(perm.begin(), perm.end()));
}

std::pair<std::vector<std::vector<double>>, double> gmm_1d::site_resp_constrained_by_site(const std::vector<std::vector<double>> &logA) const {
  const int m = logA.size();
  const int G = logA[0].size();

  if (m > G)
    throw std::runtime_error("Site has more variants than components");

  // Generate assignments of G components to m variants
  std::vector<std::vector<int>> assignments;
  generate_assignments(G, m, assignments);

  const int A = assignments.size();

  // log weights per assignment
  std::vector<double> log_w(A);

  for (int a = 0; a < A; ++a) {
    double s = 0.0;
    for (int j = 0; j < m; ++j) {
      s += logA[j][assignments[a][j]];
    }
    log_w[a] = s;
  }

  // Marginal likelihood of site data
  double log_Z = log_sum_exp(log_w);

  std::vector<double> probs(A);
  for (int a = 0; a < A; ++a) {
    probs[a] = std::exp(log_w[a] - log_Z);
  }

  // responsibilities under constraint
  std::vector<std::vector<double>> resp(m, std::vector<double>(G, 0.0));

  for (int a = 0; a < A; ++a) {
    double p = probs[a];
    for (int j = 0; j < m; ++j) {
      int g = assignments[a][j];
      resp[j][g] += p;
    }
  }

  return {resp, log_Z};
}

double gmm_1d::e_step_1d(const std::vector<double> &x, const std::vector<int> &site_id, std::vector<std::vector<double>> &resp) const {
  const int N = x.size();
  const int G = this->n_components;

  resp.assign(N, std::vector<double>(G, 0.0));

  // Compute logA
  std::vector<std::vector<double>> logA(N, std::vector<double>(G));

  for (int i = 0; i < N; ++i) {
    for (int g = 0; g < G; ++g) {
      logA[i][g] =
          std::log(std::max(this->weights[g], 1e-300)) + // avoid log(0) for weights
          log_normal_1d(x[i], this->means[g], this->vars[g]);
    }
  }

  // Group indices by site
  std::unordered_map<int, std::vector<int>> sites;
  for (int i = 0; i < N; ++i) {
    sites[site_id[i]].push_back(i);
  }

  double logL = 0.0;

  // Process each site independently
  for (const auto& kv : sites) {
    const std::vector<int>& idxs = kv.second;
    const int m = idxs.size();

    std::vector<std::vector<double>> site_logA(m, std::vector<double>(G));
    for (int j = 0; j < m; ++j) {
      site_logA[j] = logA[idxs[j]];
    }

    auto result = site_resp_constrained_by_site(site_logA);
    std::vector<std::vector<double>> site_resp = result.first;
    double site_logZ = result.second;

    for (int j = 0; j < m; ++j) {
      resp[idxs[j]] = site_resp[j];
    }

    logL += site_logZ;
  }

  return logL;
}

void gmm_1d::m_step_1d(const std::vector<double> &x, const std::vector<std::vector<double>> &resp) {
  const size_t N = resp.size();
  if (N == 0) {
    throw std::runtime_error("m_step_1d: empty resp");
  }
  const size_t G = resp[0].size();

  // L = resp.sum(axis=0)
  std::vector<double> L(G, 0.0);
  for (size_t i = 0; i < N; ++i) {
    if (resp[i].size() != G) {
      throw std::runtime_error("m_step_1d: inconsistent resp dimensions");
    }
    for (size_t g = 0; g < G; ++g) {
      L[g] += resp[i][g];
    }
  }

  this->weights.assign(G, 0.0);
  this->means.assign(G, 0.0);
  this->vars.assign(G, 0.0);

  // Precompute for resurrection of components with no weight
  std::vector<double> row_sums(N, 0.0);
  for (size_t i = 0; i < N; ++i) {
    row_sums[i] = std::accumulate(resp[i].begin(), resp[i].end(), 0.0);
  }

  // Precompute global variance
  double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / N;
  double var_x = 0.0;
  for (size_t i = 0; i < N; ++i) {
    double d = x[i] - mean_x;
    var_x += d * d;
  }
  var_x /= N;

  for (size_t g = 0; g < G; ++g) {
    if (L[g] == 0.0) {
      // Reinitialize component if likelihood is 0.
      size_t idx = std::distance(row_sums.begin(), std::max_element(row_sums.begin(), row_sums.end()));

      this->means[g] = x[idx];
      this->vars[g] = var_x;
      weights[g] = this->weight_floor;
    } else {
      weights[g] = L[g] / static_cast<double>(N);

      double mu = 0.0;
      for (size_t i = 0; i < N; ++i) {
        mu += resp[i][g] * x[i];
      }
      mu /= L[g];
      this->means[g] = mu;

      double var = 0.0;
      for (size_t i = 0; i < N; ++i) {
        double d = x[i] - mu;
        var += resp[i][g] * d * d;
      }
      var /= L[g];
      this->vars[g] = var;
    }
  }

  // set minimum var and weight
  for (size_t g = 0; g < G; ++g) {
    this->vars[g] = std::max(this->vars[g], this->var_floor);
    this->weights[g] = std::max(this->weights[g], this->weight_floor);
  }

  // renormalize weights
  double wsum = std::accumulate(weights.begin(), weights.end(), 0.0);
  for (size_t g = 0; g < G; ++g) {
    weights[g] /= wsum;
  }
}

void gmm_1d::initialize_k_means_1d(const std::vector<double> &x, int K, int n_iter, int seed) {
  const size_t N = x.size();
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> uni(0, N - 1);

  std::vector<double> centers(K);
  for (size_t k = 0; k < K; ++k) {
    centers[k] = x[uni(rng)];
  }

  std::vector<size_t> labels(N);

  for (size_t it = 0; it < n_iter; ++it) {
    // assignment
    for (size_t i = 0; i < N; ++i) {
      double best = std::numeric_limits<double>::infinity();
      size_t best_k = 0;
      for (size_t k = 0; k < K; ++k) {
        double d = x[i] - centers[k];
        double dist = d * d;
        if (dist < best) {
          best = dist;
          best_k = k;
        }
      }
      labels[i] = best_k;
    }

    // update
    std::vector<double> sum(K, 0.0);
    std::vector<size_t> count(K, 0);
    for (size_t i = 0; i < N; ++i) {
      sum[labels[i]] += x[i];
      count[labels[i]]++;
    }
    for (size_t k = 0; k < K; ++k) {
      if (count[k] > 0) {
        centers[k] = sum[k] / count[k];
      }
    }
  }

  this->means = centers;
}

bool gmm_1d::fit(const std::vector<double> &x, const std::vector<int> &sites, std::vector<double>& logL_history, int n_iter,  double tolerance, bool adaptive, unsigned int seed) {
  const size_t N = x.size();
  if (sites.size() != N) {
    throw std::runtime_error("em_gmm_1d: x and sites size mismatch");
  }
  const int G = this->n_components;
  initialize_k_means_1d(x, G, 10, seed);
  this->vars.assign(G, this->var_floor);
  this->weights.assign(G, 1.0 / static_cast<double>(G));
  std::cerr << "k-means init\nmeans: ";
  for (double v : this->means) std::cerr << v << " ";
  std::cerr << "\nvars: ";
  for (double v : this->vars) std::cerr << v << " ";
  std::cerr << "\nweights: ";
  for (double v : this->weights) std::cerr << v << " ";
  std::cerr << "\n";

  logL_history.reserve(n_iter);

  // Count unique sites
  std::vector<int> unique_sites = sites;
  std::sort(unique_sites.begin(), unique_sites.end());
  auto sites_end = std::unique(unique_sites.begin(), unique_sites.end());
  const int n_unique_sites = std::distance(unique_sites.begin(), sites_end);

  // Set convergence tolerance
  const double conv_tol = (tolerance < 0.0) ? std::numeric_limits<double>::epsilon() : tolerance;
  double old_avg_logL = -std::numeric_limits<double>::infinity();

  // EM loop
  for (size_t it = 0; it < n_iter; ++it) {

    std::vector<std::vector<double>> resp;

    // E-step
    double logL = e_step_1d(x, sites, resp);

    // M-step
    m_step_1d(x, resp);

    //TESTLINES
    if(adaptive){
      this->vars.assign(G, 0.05); 
    }
    // Average log-likelihood per site to check for convergence
    const double avg_logL = logL / static_cast<double>(n_unique_sites);

    // Check for non-finite log-likelihood
    if (std::isfinite(avg_logL) == false) {
      std::cerr << "ERROR: Non-finite log-likelihood at iteration " << it << "\n";
      return false;
    }

    // Check for convergence
    if (it > 0) {
      const double delta = std::abs(old_avg_logL - avg_logL);

      if (delta <= conv_tol) {
        std::cerr << "Converged at iteration " << it
                  << "  avg_logL=" << avg_logL
                  << "  delta=" << delta << "\n";
        logL_history.push_back(logL);
        break;
      }
    }

    old_avg_logL = avg_logL;

    logL_history.push_back(logL);

    // Logging
    std::cerr
        << "iter " << it
        << "  logL=" << logL
        << "  avg_logL=" << avg_logL;
    if (it > 0) {
      std::cerr << "  delta=" << std::abs(old_avg_logL - avg_logL);
    }
    std::cerr << "  weights=";
    for (double v : this->weights) std::cerr << v << " ";
    std::cerr << " mus=";
    for (double v : this->means) std::cerr << v << " ";
    std::cerr << " vars=";
    for (double v : this->vars) std::cerr << v << " ";
    std::cerr << "\n";
  }
  return true;
}

bool gmm_1d::predict(const std::vector<double>& x, const std::vector<int>& sites, std::vector<int>& assigned_components, std::vector<std::vector<double>>& marginal_posterior_probabilities) const {
  const int N = (int)x.size();
  const int G = this->n_components;

  if (sites.size() != N)
    throw std::runtime_error("predict(): size mismatch");

  assigned_components.resize(N);
  marginal_posterior_probabilities.resize(N);
  marginal_posterior_probabilities.assign(N, std::vector<double>(G, 0.0));

  // Calculate logA[i][g] = log pi_g + log N(x_i | g)
  std::vector<std::vector<double>> logA(N, std::vector<double>(G));

  for (int i = 0; i < N; ++i) {
    for (int g = 0; g < G; ++g) {
      logA[i][g] =
          std::log(std::max(weights[g], 1e-300)) +
          log_normal_1d(x[i], means[g], vars[g]);
    }
  }

  // Group variants by site
  std::unordered_map<int, std::vector<int>> sites_variants;
  for (int i = 0; i < N; ++i)
    sites_variants[sites[i]].push_back(i);

  for (const auto& kv : sites_variants) {
    const std::vector<int>& idxs = kv.second;
    const int m = (int)idxs.size();

    if (m > G)
      throw std::runtime_error("predict(): site size > #components");

    // site_logA[j][g]
    std::vector<std::vector<double>> site_logA(
        m, std::vector<double>(G));

    for (int j = 0; j < m; ++j)
      site_logA[j] = logA[idxs[j]];

    std::vector<std::vector<int>> assignments;
    generate_assignments(G, m, assignments);

    const int K = (int)assignments.size();
    std::vector<double> loglik(K);

    // compute log-likelihood for each assignment
    for (int k = 0; k < K; ++k) {
      double s = 0.0;
      for (int j = 0; j < m; ++j)
        s += site_logA[j][assignments[k][j]];
      loglik[k] = s;
    }

    const double Z = log_sum_exp(loglik);

    // Calculate marginal posterior probabilities
    for (int k = 0; k < K; ++k) {
      const double w = std::exp(loglik[k] - Z);
      for (int j = 0; j < m; ++j) {
        int i = idxs[j];
        int g = assignments[k][j];
        marginal_posterior_probabilities[i][g] += w;
      }
    }

    int k_star = std::max_element(loglik.begin(), loglik.end()) - loglik.begin();

    for (int j = 0; j < m; ++j)
      assigned_components[idxs[j]] = assignments[k_star][j];
  }
  return true;
}

double gmm_1d::get_log_likelihood(const std::vector<double> &x, const std::vector<int> &sites) const {
  if (x.size() != sites.size())
    throw std::runtime_error("log_likelihood(): size mismatch");

  // Reusing e_step_1d with redundant resp
  std::vector<std::vector<double>> resp;
  const double ll = this->e_step_1d(x, sites, resp);

  if (!std::isfinite(ll))
    throw std::runtime_error("log_likelihood(): non-finite log-likelihood");

  return ll;
}

double gmm_1d::get_bic(const std::vector<double> &x, const std::vector<int> &sites) const {
  if (x.size() != sites.size())
    throw std::runtime_error("bic(): size mismatch");

  const int G = this->n_components;
  if (G <= 0) throw std::runtime_error("bic(): invalid n_components");

  // Calcualte N = number of independent sites
  std::vector<int> unique_sites = sites;     // copy
  std::sort(unique_sites.begin(), unique_sites.end());
  auto it = std::unique(unique_sites.begin(), unique_sites.end());
  const int n_unique_sites = std::distance(unique_sites.begin(), it);

  const double n = static_cast<double>(n_unique_sites);
  if (n <= 1)
    throw std::runtime_error("bic(): need at least 2 sites for BIC");

  const double k = static_cast<double>(3 * G - 1);

  const double logL = this->get_log_likelihood(x, sites);

  // BIC = -2*logL + k*log(n)
  return -2.0 * logL + k * std::log(n);
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
int gmm_1d::get_distinct_components_count(const std::vector<int>& sites, double min_bd_threshold) {
  const int G = n_components;
  if (G <= 0) return 0;

  // Get minimum number of components from site constraint
  int m_max = 0;
  {
    std::vector<int> s = sites;
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

  std::vector<int> active;
  active.reserve(G);
  for (int g = 0; g < G; ++g) {
    if (weights[g] > DEFAULT_WEIGHT_FLOOR)
      active.push_back(g);
  }

  if (active.empty())
    return m_max;

  std::vector<std::vector<int>> clusters;
  clusters.reserve(active.size());
  for (int g : active)
    clusters.push_back({g});

  // Merge clusters until all BD > threshold for all clusters
  bool merged;
  do {
    merged = false;

    // Compare clusters pairwise
    for (int i = 0; i < clusters.size() && !merged; ++i) {
      for (int j = i + 1; j < clusters.size() && !merged; ++j) {

        bool can_merge = true;

        for (int gi : clusters[i]) {
          for (int gj : clusters[j]) {
            double d = calculate_bhattacharyya_distance_1d(means[gi], vars[gi], means[gj], vars[gj]);
            if (d > min_bd_threshold) {
              can_merge = false;
              break;
            }
          }
          if (!can_merge)
            break;
        }

        if (can_merge) {
          // Merge clusters
          clusters[i].insert(
              clusters[i].end(),
              clusters[j].begin(),
              clusters[j].end()
          );
          clusters.erase(clusters.begin() + j);
          merged = true;
        }
      }
    }
  } while (merged);

  std::vector<double> m_sigmoid(means.size());
  sigmoid_transform(means, m_sigmoid);

  for (int c = 0; c < clusters.size(); ++c) {

    int rep = clusters[c][0];
    double max_w = weights[rep];

    for (int g : clusters[c]) {
      if (weights[g] > max_w) {
        max_w = weights[g];
        rep = g;
      }
    }
    merged_means.push_back(m_sigmoid[rep]);
    merged_vars.push_back(vars[rep]);
    merged_weights.push_back(weights[rep]);

    std::cerr << "distinct cluster " << c
              << ": mean=" << means[rep]
              << " mean in frequency space=" << m_sigmoid[rep]
              << " var=" << vars[rep]
              << " weight=" << weights[rep]
              << "\n";
  }

  int g_unique = static_cast<int>(clusters.size());
  if (g_unique < m_max)
    g_unique = m_max;

  return g_unique;
}
