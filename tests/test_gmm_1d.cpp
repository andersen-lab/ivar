#include "../src/gmm_1d.h"
#include "../src/include/armadillo"
#include <random>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <sstream>

void read_in_simulated_data(std::string fname, std::vector<double>& frequencies, std::vector<int>& sites, std::vector<uint32_t> &depths, float invariant_threshold) {
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

    if(alt_freq <= invariant_threshold && alt_freq >= (1 - invariant_threshold) && total_dp >= 10 && alt_qual >= 20) {
      frequencies.push_back(alt_freq);
      sites.push_back(std::stoi(row.at(1)));
      depths.push_back(std::stoi(row.at(7)));
    }

  }
}

int main(int argc, char* argv[]) {

  std::vector<double> x;
  std::vector<int> sites;
  std::vector<uint32_t> depths;
  std::string input = argv[1];
  std::string prefix = argv[2];
  double eps = 1e-6;

  std::vector<double> y;

  sites.clear();
  read_in_simulated_data(input, x, sites, depths, 0.97);

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

  // Logit transform
//  std::vector<double> x_logit;
//  gmm_1d::logit_transform(x, x_logit, eps);

  std::ofstream out(prefix + "_gmm_1d_results.txt");

  out << "Replicate\tComponents\tBIC\tDistinct_Components\tMeans\tVariances\tWeights\n";
  // Fit GMM
  for(int b = 0; b < 100; b++) {
    for(int i = 4; i < 8; i++)  {
      std::vector<double> logL_history;
      gmm_1d model(i);
//      model.set_use_half_normal_for_noise(false);
      model.fit(
          x,
          sites,
          logL_history,
          depths,
          20,
          1e-6
      );

      std::cerr << "Distinct components: " << i << " " << model.get_distinct_components_count(sites) << "\n";

//      std::vector<double> m_sigmoid;
      std::vector<double> means = model.get_means();
      std::vector<double> vars = model.get_vars();
      std::vector<double> weights = model.get_weights();

//      gmm_1d::sigmoid_transform(means, m_sigmoid, eps);

      out << b << "\t"
          << i << "\t"
          << model.get_bic(x, sites) << "\t"
          << model.get_distinct_components_count(sites) << "\t";

      for(int j = 0; j < means.size(); j++){
        out << means[j];
        if (j < means.size() - 1){
          out << ",";
        }
      }
      out << "\t";

      for(int j = 0; j < vars.size(); j++) {
        out << vars[j];
        if (j < vars.size() - 1) {
          out << ",";
        }
      }
      out << "\t";

      for(int j = 0; j < weights.size(); j++) {
        out << weights[j];
        if (j < weights.size() - 1) {
          out << ",";
        }
      }
      out << "\n";
//
//      for(int j = 0; j < m_sigmoid.size(); j++) {
//        out << m_sigmoid[j];
//        if (j < m_sigmoid.size() - 1) {
//          out << ",";
//        }
//      }
//      out << "\n";

    }
  }

  out.close();

  return 0;
}

