#include "../src/site_state_aggregator_stats.h"

bool check(const amplicon_variation_data &av, const std::vector<double> &expected_frequencies,
           const std::vector<uint32_t> &expected_weights) {
  if(av.num_amplicons != (int)expected_frequencies.size()) {
    std::cerr << "Expected " << expected_frequencies.size() << " amplicons, got " << av.num_amplicons << std::endl;
    return false;
  }
  if(av.frequencies.size() != expected_frequencies.size() || av.weights.size() != expected_weights.size()) {
    std::cerr << "Frequency/weight vectors are not aligned with the amplicon list" << std::endl;
    return false;
  }
  for(uint32_t i = 0; i < expected_frequencies.size(); i++){
    if(av.frequencies[i] != expected_frequencies[i]){
      std::cerr << "Frequencies do not match expected values" << std::endl;
      return false;
    }
    if(av.weights[i] != expected_weights[i]){
      std::cerr << "Weights do not match expected values" << std::endl;
      return false;
    }
  }
  return true;
}

int main() {
  ITNode *amp1 = new ITNode(Interval(0, 25));
  ITNode *amp2 = new ITNode(Interval(15, 40));

  // Test with the state present on two amplicons
  site_state_aggregator_stats ssas("A");

  for(int i = 0; i < 10; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp1);
    ssas.add_site(ss);
  }

  for(int i = 0; i < 5; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp2);
    ssas.add_site(ss);
  }

  amplicon_variation_data av = ssas.calculate_amplicon_variation({ {amp1, 10}, {amp2, 10} });
  if(!check(av, {1.0, 0.5}, {10, 10})) return -1;

  // Test with the state present on one of two covering amplicons. The amplicon it is
  // absent from must still be reported, at frequency 0, so dropout stays visible
  site_state_aggregator_stats ssas2("A");

  for(int i = 0; i < 15; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp1);
    ssas2.add_site(ss);
  }

  av = ssas2.calculate_amplicon_variation({ {amp1, 20}, {amp2, 10} });
  if(!check(av, {0.75, 0.0}, {20, 10})) return -1;

  // Reads that could not be assigned to a single amplicon are keyed under nullptr and
  // must be excluded rather than reported as an amplicon
  site_state_aggregator_stats ssas3("A");

  for(int i = 0; i < 10; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(amp1);
    ssas3.add_site(ss);
  }

  for(int i = 0; i < 15; i++){
    site_state ss("A", 20, 20, NUCLEOTIDE);
    ss.set_amplicon(nullptr);
    ssas3.add_site(ss);
  }

  av = ssas3.calculate_amplicon_variation({ {amp1, 20}, {amp2, 10}, {nullptr, 15} });
  if(!check(av, {0.5, 0.0}, {20, 10})) return -1;

  // Amplicons with no depth at the position are not reported
  av = ssas3.calculate_amplicon_variation({ {amp1, 20}, {amp2, 0} });
  if(!check(av, {0.5}, {20})) return -1;

  return 0;
}
