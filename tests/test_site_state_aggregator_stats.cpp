#include "../src/site_state_aggregator_stats.h"

int main() {
  std::vector<site_state> site_states;

  ITNode *amp1 = new ITNode(Interval(0, 25));
  ITNode *amp2 = new ITNode(Interval(15, 40));

  // Test with two amplicons
  site_state_aggregator_stats ssas("A");

  for(int i = 0; i < 10; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp1);
    ssas.add_site(ss);
  }

  for(int i = 0; i < 5; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp2);
    site_states.emplace_back(ss);
    ssas.add_site(ss);
  }

  amplicon_variation_data av = ssas.calculate_amplicon_variation({ {amp1, 10}, {amp2, 10} });

  if(av.num_amplicons != 2 || av.stdev != 0.25 || !av.amplicon_frequency_variation) {
    std::cerr << "Expected 2 amplicons with stdev 0.25 and variation true" << std::endl;
    return -1;
  }

  std::vector<double> expected_frequencies = {1.0, 0.5};
  for(int i = 0; i < av.frequencies.size(); i++){
    if(av.frequencies[i] != expected_frequencies[i]){
      std::cerr << "Frequencies do not match expected values" << std::endl;
      return -1;
    }
  }
  std::vector<double> expected_weights = {10, 10};
  for(int i = 0; i < av.weights.size(); i++){
    if(av.weights[i] != expected_weights[i]){
      std::cerr << "Weights do not match expected values" << std::endl;
      return -1;
    }
  }

  // Test with one amplicon
  site_state_aggregator_stats ssas2("A");
  expected_frequencies.clear();
  expected_weights.clear();

  for(int i = 0; i < 10; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp1);
    ssas2.add_site(ss);
  }

  for(int i = 0; i < 5; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp1);
    site_states.emplace_back(ss);
    ssas2.add_site(ss);
  }

  av = ssas2.calculate_amplicon_variation({ {amp1, 20}, {amp2, 10} });

  if(av.num_amplicons != 1 || av.stdev != 0 || av.amplicon_frequency_variation || av.frequencies.size() != 0 || av.weights.size() != 0) {
    std::cerr << "Expected 1 amplicons with stdev 0 and variation false" << std::endl;
    return -1;
  }

  // Test with null
  site_state_aggregator_stats ssas3("A");
  expected_frequencies.clear();
  expected_weights.clear();

  for(int i = 0; i < 10; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(nullptr);
    ssas2.add_site(ss);
  }

  for(int i = 0; i < 5; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(nullptr);
    site_states.emplace_back(ss);
    ssas2.add_site(ss);
  }

  av = ssas2.calculate_amplicon_variation({ {amp1, 20}, {amp2, 10} });

  if(av.num_amplicons != 1 || av.stdev != 0 || av.amplicon_frequency_variation || av.frequencies.size() != 0 || av.weights.size() != 0) {
    std::cerr << "Expected 1 amplicons with stdev 0 and variation false" << std::endl;
    return -1;
  }

  return 0;
}