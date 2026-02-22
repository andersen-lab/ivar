#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <unordered_map>

#include "../src/bootstrap_variants.h"
#include "../src/gmm_1d.h"

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
  std::vector<uint32_t> depths;
  std::vector<uint32_t> total_depths;
  std::string input = argv[1];
  std::string prefix = argv[2];
  double eps = 1e-6;

  std::vector<double> y;

  std::vector<uint32_t> sites;
  read_in_simulated_data(input, x, sites, depths, total_depths, 0.97);

  std::ofstream out(prefix + "_gmm_1d_results.txt");

  out << "Replicate\tComponents\tBIC\tDistinct_Components\tMeans\tVariances\tWeights\n";
  // Fit GMM
  for(int b = 0; b < 1; b++) {
//    bootstrap_variants bootstrap(sites, depths, total_depths);

//    std::vector<uint32_t> sampled_sites;
//    std::vector<uint32_t> sampled_depths;
//    std::vector<double> sampled_frequencies;

//    bootstrap.sample(sampled_sites, sampled_depths, sampled_frequencies, counts.size());

    gmm_1d model(12, 42);  // seed matches bootstrap replicate

    model.fit(x);
    std::vector<int> labels = model.predict(x);

    std::vector<int> component_indices = model.get_effective_components(labels);
    std::cerr << "VB effective components: " << component_indices.size() << "\n";

    std::cerr << "VB effective means: ";
    std::vector<double> eff_means = model.get_effective_means(component_indices);
    for(int i = 0; i < eff_means.size(); i++) {
      std::cerr << eff_means[i] << " ";
    }
    std::cerr << "\n";

    std::cerr << "VB effective vars: ";
    std::vector<double> eff_vars = model.get_effective_vars(component_indices);
    for(int i = 0; i < eff_vars.size(); i++) {
      std::cerr << eff_vars[i] << " ";
    }
    std::cerr << "\n";

    std::cerr << "VB effective weights: ";
    std::vector<double> eff_weights = model.get_effective_weights(component_indices);
    for(int i = 0; i < eff_weights.size(); i++) {
      std::cerr << eff_weights[i] << " ";
    }
    std::cerr << "\n";

    std::cerr << "ELBO history: ";
    for(int i = 0; i < model.get_elbo_history().size(); i++) {
      std::cerr << model.get_elbo_history()[i] << ", ";
    }

    std::vector<double> means   = model.get_means();
    std::vector<double> vars    = model.get_variances();
    std::vector<double> weights = model.get_weights();

    out << b << "\t"
        << "VB\t"
        << "NA\t"  // no BIC for VB
        << component_indices.size() << "\t";

    for (size_t j = 0; j < means.size(); j++) {
      out << means[j];
      if (j < means.size() - 1) out << ",";
    }
    out << "\t";

    for (size_t j = 0; j < vars.size(); j++) {
      out << vars[j];
      if (j < vars.size() - 1) out << ",";
    }
    out << "\t";

    for (size_t j = 0; j < weights.size(); j++) {
      out << weights[j];
      if (j < weights.size() - 1) out << ",";
    }
    out << "\n";
  }

  out.close();

  return 0;
}

