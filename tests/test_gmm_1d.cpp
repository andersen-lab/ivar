#include "../src/gmm_1d.h"
#include "../src/bootstrap_variants.h"
#include "../src/solve_clustering.h"
#include "../src/genomic_position.h"
#include <random>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <sstream>

void read_in_simulated_data(std::string fname, std::vector<double>& frequencies, std::vector<uint32_t>& sites, std::vector<uint32_t> &depths, std::vector<uint32_t> &total_depths, std::vector<std::string> &nucs, float invariant_threshold) {
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
      nucs.push_back(row.at(3));
      depths.push_back(std::stoi(row.at(7))); // Using ALT_DP for now
      total_depths.push_back(total_dp);
    }

  }
}

int get_max_variants_per_site(const std::vector<uint32_t> &sites, const std::vector<std::string> nucs) {
  std::unordered_map<uint32_t, int> counts;
  for (size_t i = 0; i < sites.size(); i++) {
    bool contains = nucs[i].find("+") != std::string::npos;
    if (contains) {
      continue;
    }
    //contains = nucs[i].find("-") != std::string::npos; 
    counts[sites[i]]++;
    if(counts[sites[i]] > 2) {
      std::cerr << sites[i] << std::endl;
    }
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
  std::vector<std::string> nucs;
  std::string input = argv[1];
  std::string prefix = argv[2];
  double eps = 1e-6;

  std::vector<double> y;

  sites.clear();
  read_in_simulated_data(input, x, sites, depths, total_depths, nucs, 0.97);
  int N_min = get_max_variants_per_site(sites, nucs);

  x.clear();
  sites.clear();
  depths.clear();
  total_depths.clear();
  nucs.clear(); 
  read_in_simulated_data(input, x, sites, depths, total_depths, nucs, 0.99);

  std::cerr << "N_min: " << N_min << std::endl; 
  std::ofstream out(prefix + "_gmm_1d_results.txt");

  out << "N_min\tOriginal_Components\tBIC\tDistinct_Components\tMeans\tVariances\tWeights\tPass_Subset\n";
  uint32_t largest_passing_n = std::numeric_limits<uint32_t>::min();
 
  for(uint32_t i=N_min+2; i < N_min+2+4; i++){
    std::cerr << "\n"; 
    std::vector<double> logL_history;
    gmm_1d model(i, 0);
    model.set_use_half_normal_for_noise(true);
    model.fit(
        x,
        sites,
        logL_history,
        depths,
        20,
        1e-6
    );

    std::cerr << "Original components: " << i << " Distinct components: " << model.get_distinct_components_count(sites) << "\n";
    std::vector<double> logL_history2;
    gmm_1d model2(model.get_distinct_components_count(sites), 0);
    model2.set_use_half_normal_for_noise(true);
    model2.fit(
        x,
        sites,
        logL_history,
        depths,
        20,
        1e-6
    );

    std::vector<double> means = model2.get_means();
    std::vector<double> vars = model2.get_vars();
    std::vector<double> weights = model.get_weights();

    std::vector<std::vector<double>> solution_sets;
    bool pass_subset = subset_sum(means, solution_sets, 0.10); 
    if(pass_subset){
      if(means.size() > largest_passing_n){
        largest_passing_n = means.size();
      }
    }
    std::cerr << "pass subset " << pass_subset << "\n";

    out << N_min << "\t"
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
    if(pass_subset){
      out << "\t" << "Yes";
    } else {
      out << "\t" << "No";
    }
    out << "\n";
  }
  out.close();
  std::cerr << "largest passing number of components: " << largest_passing_n << "\n";
  return 0;
}

