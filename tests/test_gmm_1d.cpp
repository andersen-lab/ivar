#include "../src/gmm_1d.h"
#include "../src/bootstrap_variants.h"
#include <random>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <sstream>

void read_in_simulated_data(std::string fname, std::vector<double>& frequencies, std::vector<uint32_t>& sites, std::vector<uint32_t> &depths, std::vector<uint32_t> &total_depths, float invariant_threshold) {
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

    uint32_t total_dp = std::stoi(row.at(11));
    uint8_t alt_qual = std::stoi(row.at(9));

    if(alt_freq <= invariant_threshold && alt_freq >= (1 - invariant_threshold) && total_dp >= 10 && alt_qual >= 20) {
      frequencies.push_back(alt_freq);
      sites.push_back(std::stoi(row.at(1)));
//      depths.push_back(std::stoi(row.at(21)));
      depths.push_back(std::stoi(row.at(7))); // Using ALT_DP for now
      total_depths.push_back(total_dp);
    }

  }
}

int get_max_variants_per_site(const std::vector<uint32_t> &sites) {
  std::unordered_map<uint32_t, int> counts;
  for (uint32_t v : sites) {
    counts[v]++;
  }
  int N_min = -std::numeric_limits<double>::infinity();
  for (const auto& [key, value] : counts) {
    N_min = std::max(N_min, value);
  }
  return N_min;
}

int main(int argc, char* argv[]) {

  std::vector<double> x;
  std::vector<uint32_t> sites;
  std::vector<uint32_t> depths;
  std::vector<uint32_t> total_depths;
  std::string input = argv[1];
  std::string prefix = argv[2];
  double eps = 1e-6;

  std::vector<double> y;

  sites.clear();
  read_in_simulated_data(input, x, sites, depths, total_depths, 0.97);


  int N_min = get_max_variants_per_site(sites);

  // Hacky way to get unique number of sites
  std::unordered_map<uint32_t, int> counts;
  for (uint32_t v : sites) {
    counts[v]++;
  }

  // Logit transform
//  std::vector<double> x_logit;
//  gmm_1d::logit_transform(x, x_logit, eps);

  std::ofstream out(prefix + "_gmm_1d_results.txt");

  out << "Replicate\tComponents\tBIC\tDistinct_Components\tMeans\tVariances\tWeights\n";
  // Fit GMM
  for(int b = 0; b < 100; b++) {
    bootstrap_variants bootstrap(sites, depths, total_depths);

    std::vector<uint32_t> sampled_sites;
    std::vector<uint32_t> sampled_depths;
    std::vector<double> sampled_frequencies;

    bootstrap.sample(sampled_sites, sampled_depths, sampled_frequencies, counts.size());

    std::cerr << "Sampled sites size: " << sampled_sites.size() << std::endl;
    std::cerr << "Sampled depths size: " << sampled_depths.size() << std::endl;

    N_min = get_max_variants_per_site(sampled_sites);
    N_min = std::max(N_min, 2);

    for(int i = N_min + 2; i < N_min + 6; i++)  {

      std::vector<double> logL_history;
      gmm_1d model(i, b);
      std::cerr << "Seed: " << model.get_seed() << std::endl;
      std::cerr << "Number of components: " << N_min << std::endl;

//      model.set_use_half_normal_for_noise(false);
      model.fit(
          sampled_frequencies,
          sampled_sites,
          logL_history,
          sampled_depths,
          20,
          1e-6
      );

      std::cerr << "Distinct components: " << i << " " << model.get_distinct_components_count(sampled_sites) << "\n";

//      std::vector<double> m_sigmoid;
      std::vector<double> means = model.get_means();
      std::vector<double> vars = model.get_vars();
      std::vector<double> weights = model.get_weights();

//      gmm_1d::sigmoid_transform(means, m_sigmoid, eps);

      out << b << "\t"
          << i << "\t"
          << model.get_bic(sampled_frequencies, sampled_sites) << "\t"
          << model.get_distinct_components_count(sampled_sites) << "\t";

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

