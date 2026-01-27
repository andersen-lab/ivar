#include "../src/bootstrap_variants.h"
#include <iostream>

int main() {
  std::vector<uint32_t> sites = {1, 1, 2, 2, 2, 3, 3};
  std::vector<uint32_t> depths = {5, 5, 25, 70, 5, 5, 45};
  std::vector<uint32_t> total_depths= {10, 10, 100, 100, 100, 50, 50};

  bootstrap_variants bootstrap(sites, depths, total_depths);

  std::vector<uint32_t> sampled_variants;
  std::vector<uint32_t> sampled_depths;
  std::vector<double> sampled_frequencies;
  bootstrap.sample(sampled_variants, sampled_depths, sampled_frequencies, 3);

  for(int i = 0; i < sampled_variants.size(); i++) {
      std::cout << "Sampled Variant Index: " << sampled_variants[i]
                << ", Depth: " << sampled_depths[i]
                << ", Frequency: " << sampled_frequencies[i] << std::endl;
  }

  return 0;
}
