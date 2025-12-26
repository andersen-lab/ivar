#include "../src/gmm_1d.h"
#include <random>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

double log_sum_exp_naive(const double* x, int n) {
  double sum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sum += std::exp(x[i]);
  }
  return std::log(sum);
}

void read_in_variants(std::string fname, std::vector<double>& frequencies, std::vector<int>& sites) {
  std::ifstream infile(fname);
  std::string line;

  while (std::getline(infile, line)) {
    std::istringstream ss(line);
    std::string cell;
    std::vector<std::string> row;

    while (std::getline(ss, cell, '\t')) {
      row.push_back(cell);
    }

    if(row.at(10).compare("ALT_FREQ") == 0)
      continue;
    double alt_freq = std::stod(row.at(10));
    int total_dp = std::stoi(row.at(11));

    float invariant_threshold = 0.95;

    if(alt_freq < invariant_threshold && alt_freq > (1 - invariant_threshold) && total_dp >= 10) {
      frequencies.push_back(alt_freq);
      sites.push_back(std::stoi(row.at(1)));
    }

  }
}

int main() {

//  double x[] = {-107498545664.81509, -9363075829.8171787};
//  double y = gmm_1d::log_sum_exp(x, 2);
//  double y_expected = log_sum_exp_naive(x, 2);

//  std::cerr << y << " " <<  y_expected << std::endl;

//  const std::vector<double> means = {0.2, 0.4, 0.6, 0.8};
//  const std::vector<double> sds = {0.05, 0.1, 0.1, 0.05};
//  const std::vector<int> points_per_cluster = {20, 10, 10, 10};
//  std::vector<double> x;
//  std::vector<int> sites;
//  std::mt19937 rng(112358);
//
//  const int total_points = std::accumulate(points_per_cluster.begin(), points_per_cluster.end(), 0);
//  x.resize(total_points);
//  sites.resize(total_points);
//
//  int offset = 0;
//  for (int k = 0; k < means.size(); ++k) {
//    std::normal_distribution<double> dist(means[k], sds[k]);
//
//    for (int i = 0; i < points_per_cluster[k]; ++i) {
//      x[offset] = dist(rng);
//      sites[offset] = i;
//      ++offset;
//    }
//  }
//
//  // Clip to (0,1) to avoid logit issues
  const double eps = 1e-6;
//  for (double& v : x) {
//    if (v <= eps) v = eps;
//    if (v >= 1.0 - eps) v = 1.0 - eps;
//  }
//
//  std::ofstream out("/Users/karthik/Documents/code/saga/output_x.tsv");
//  out << "data" << '\n';
//  for (double v : x) {
//    out << v << '\n';
//  }
//  out.close();

  std::vector<double> x;
  std::vector<int> sites;
  int N = 4;

//  read_in_variants("/Users/karthik/tmp/ivar_chronic_paper/lineages_three_populations/var/triple_mix_2.txt", x, sites);
  read_in_variants("/Users/karthik/tmp/ivar_chronic_paper/variants/SRR23446131_var.txt", x, sites);


  // Logit transform
  std::vector<double> x_logit;
  gmm_1d::logit_transform(x, x_logit, eps);


  // Fit GMM
  std::vector<double> logL_history;
  gmm_1d model(N);
  model.fit(
      x_logit,
      sites,
      logL_history
  );

  // Print results (logit space)
  const auto& w = model.get_weights();
  const auto& m = model.get_means();
  const auto& v = model.get_vars();

  std::cout << "\nFitted GMM:\n";

  std::vector<double> m_sigmoid;
  gmm_1d::sigmoid_transform(m, m_sigmoid, eps);

  for (size_t g = 0; g < N; ++g) {
    std::cout
        << "Component " << g
        << "  mean=" << m[g]
        << "  mean in frequency space=" << m_sigmoid[g]
        << "  weight=" << w[g]
        << "  var=" << v[g]
        << "\n";
  }

  std::cout << "Unique components: " << model.get_distinct_components_counts(sites) << "\n";

  std::vector<int> assigned_components;
  std::vector<std::vector<double>> marginal_posterior_probabilities;
  model.predict(
      x_logit,
      sites,
      assigned_components,
      marginal_posterior_probabilities
  );
  std::cerr << "BIC: "<< model.get_bic(x_logit, sites) << "\n";

  std::ofstream out("/Users/karthik/Documents/code/saga/output_x.tsv");
  out << "data" << "\t" << "site" << "\t";
  for (int i = 0; i < N; ++i) {
    out << "comp_" << i;
    if (i < N - 1) {
      out << "\t";
    }
  }
  out << "\n";
  for (int i = 0; i < x.size(); i++) {
    out << x[i] << "\t";
    out << sites[i] << "\t";
    for (int j = 0; j < N; j++) {
      out << marginal_posterior_probabilities[i][j];
      if (j < N - 1) {
        out << "\t";
      }
    }
    out << "\n";
  }
  out.close();

  return 0;
}

