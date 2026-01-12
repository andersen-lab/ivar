#include "../src/gmm_1d.h"
#include "../src/include/armadillo"
#include <random>
#include <iostream>
#include <unordered_map>
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

void read_in_simulated_data(std::string fname, std::vector<double>& frequencies, std::vector<uint32_t>& sites, float invariant_threshold) {
  std::ifstream infile(fname);
  std::string line;
  uint32_t i=0;
  while (std::getline(infile, line)) {
    if(i ==0){
      i++;
      continue;
    }
    std::istringstream ss(line);
    std::string cell;
    std::vector<std::string> row;

    while (std::getline(ss, cell, '\t')) {
      row.push_back(cell);
    }

    double alt_freq = std::stod(row.at(20));

    double total_dp = std::stod(row.at(11));
    double alt_qual = std::stod(row.at(9));

    if(alt_freq < invariant_threshold && alt_freq > (1 - invariant_threshold) && total_dp >= 10 && alt_qual >= 20) {
      frequencies.push_back(alt_freq);
      sites.push_back(std::stoi(row.at(1)));
    }

  }
}
void read_in_universal_data(std::string fname, std::vector<double>& frequencies, std::vector<uint32_t> &sites) {
  std::ifstream infile(fname);
  std::string line;
  uint32_t i=0;
  while (std::getline(infile, line)) {
    if(i ==0){
      i++;
      continue;
    }
    std::istringstream ss(line);
    std::string cell;
    std::vector<std::string> row;

    while (std::getline(ss, cell, '\t')) {
      row.push_back(cell);
    }

    double alt_freq = std::stod(row.at(20));

    float lower_threshold = 0.50;
    float upper_threshold = 0.99;
    double total_dp = std::stod(row.at(11));
    double alt_qual = std::stod(row.at(9));

    if(alt_freq < upper_threshold && alt_freq > lower_threshold && total_dp >= 10 && alt_qual >= 20) {
      frequencies.push_back(alt_freq);
      sites.push_back(std::stoi(row.at(1)));
    }

  }
}

std::vector<double> z_score(std::vector<double> data, double mean) {
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double var = sq_sum / data.size() - (mean * mean);
    var = std::max(0.0, var);
    if(var == 0){
      var = 1e-6;
    }
    double stddev = std::sqrt(var);

    std::vector<double> z_scores;
    for (double x : data)
        z_scores.push_back((x - mean) / stddev);
    return z_scores;
}

double find_adaptive_threshold(const std::vector<double>& frequencies, const double eps, std::vector<uint32_t> sites){
  std::vector<double> x_logit;
  gmm_1d::logit_transform(frequencies, x_logit, eps);
  uint32_t N = 5;

  // Fit GMM
  std::vector<double> logL_history;
  gmm_1d model(N);
  model.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6,
      true
  );

  // Print results (logit space)
  const auto& w = model.get_weights();
  const auto& m = model.get_means();
  const auto& v = model.get_vars();

  std::vector<double> m_sigmoid;
  gmm_1d::sigmoid_transform(m, m_sigmoid, eps);
  model.get_distinct_components_count(sites);

  std::cerr << "Merged components: " << model.merged_means.size() << "\n";
  
  //extra check to make sure that the model means aren't super close together...
  uint32_t final_N = model.merged_means.size();
  gmm_1d model2(final_N);
  model2.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6,
      true
  );
  const auto& m2 = model2.get_means();
  std::vector<double> m_sigmoid2;
  gmm_1d::sigmoid_transform(m2, m_sigmoid2, eps);
  double largest_mean = 0.0;
  uint32_t largest_idx = 0;
  for(uint32_t m=0; m < m_sigmoid2.size(); ++m){
    std::cerr << "mean " << m_sigmoid2[m] << "\n";
    if(m_sigmoid2[m] > largest_mean){
      largest_mean = m_sigmoid2[m];
      largest_idx = m;
    }
  }
  std::vector<uint32_t> assigned_components;
  std::vector<std::vector<double>> marginal_posterior_probabilities;
  model2.predict(
      x_logit,
      sites,
      assigned_components,
      marginal_posterior_probabilities
  );

  std::vector<double> largest_cluster_snps;
  for(uint32_t i =0; i < frequencies.size(); i++){
    if(assigned_components[i] == largest_idx){
      if(marginal_posterior_probabilities[i][largest_idx] < 0.50){
        std::cerr << "freq " << frequencies[i] << " " << marginal_posterior_probabilities[i][largest_idx] << "\n";
      }
      largest_cluster_snps.push_back(frequencies[i]);
    }
  }

  double smallest_cluster_snp = 1.0;
  for(uint32_t k=0; k < largest_cluster_snps.size(); k++){
    if(largest_cluster_snps[k] < smallest_cluster_snp){
      smallest_cluster_snp = largest_cluster_snps[k];
    }
  }

  return(smallest_cluster_snp);
}

int main(int argc, char* argv[]) {
  const double eps = 1e-6;

  std::vector<double> x;
  std::vector<uint32_t> sites;
  std::string input = argv[1];
  std::string prefix = argv[2];

  //determine the adaptive threshold
  std::vector<double> y; 
  read_in_universal_data(input, y, sites);
  double largest_mean = find_adaptive_threshold(y, eps, sites);
  std::cerr << "Adaptive threshold: " << largest_mean << "\n";
  exit(0);
  sites.clear();
  read_in_simulated_data(input, x, sites, largest_mean);

  int n_min = 0;
  std::unordered_map<int, int> counts;
  for (int v : sites) {
      ++counts[v];
  }

  std::unordered_map<int, int>::const_iterator it =
      std::max_element(
          counts.begin(), counts.end(),
          [](const std::pair<const int, int>& a,
            const std::pair<const int, int>& b) {
              return a.second < b.second;
          }
      );

  if (it != counts.end()) {
      n_min = it->second;
  }
  int N = 0;
  std::cerr << "min n " << n_min << "\n";
  if (n_min == 3){
    N = 7;
  } else if (n_min == 4) {
    N = 9;
  } else if (n_min == 2){
    N = 5;
  } else {
    std::cerr << "Insufficient data to fit GMM\n";
    exit(0);
  }
  std::cerr << "Fitting GMM with " << N << " components\n";
  //exit(0);

  // Logit transform
  std::vector<double> x_logit;
  gmm_1d::logit_transform(x, x_logit, eps);

  // Fit GMM
  std::vector<double> logL_history;
  gmm_1d model(N);
  model.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6
  );

  // Print results (logit space)
  const auto& w = model.get_weights();
  const auto& m = model.get_means();
  const auto& v = model.get_vars();

  std::vector<double> m_sigmoid;
  gmm_1d::sigmoid_transform(m, m_sigmoid, eps);
  model.get_distinct_components_count(sites);
  std::cerr << "Merged components: " << model.merged_means.size() << "\n";

  //logging for number of components versus number of merged components
  std::ofstream merge_out(prefix + "_merged.txt");
  merge_out << "original_N\tmerged_N\n";
  merge_out << std::to_string(N) << "\t";
  merge_out << std::to_string(model.merged_means.size()) << "\n";
  merge_out.close();
  exit(0);

  int final_n = model.merged_means.size();

  gmm_1d model2(final_n);
  model2.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6
  );
  const auto& w2 = model2.get_weights();
  const auto& m2 = model2.get_means();
  const auto& v2 = model2.get_vars();
  std::vector<double> m_sigmoid2;
  gmm_1d::sigmoid_transform(m2, m_sigmoid2, eps); 

  for(auto x : m_sigmoid2){
    std::cerr << "mean second model " << x << "\n";
  }

  std::ofstream mean_out(prefix + ".txt");
  mean_out << "mean\tweight\tvar\n";
  for(uint32_t i=0; i < model.merged_means.size(); i++){
    mean_out << std::to_string(m_sigmoid2[i]) << "\t";
    mean_out << std::to_string(w2[i]) << "\t";
    mean_out << std::to_string(v2[i]) << "\n";
  }
  mean_out.close();
  exit(0);

  std::vector<uint32_t> assigned_components;
  std::vector<std::vector<double>> marginal_posterior_probabilities;
  model.predict(
      x_logit,
      sites,
      assigned_components,
      marginal_posterior_probabilities
  );

  std::ofstream out("./output_x.tsv");
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

