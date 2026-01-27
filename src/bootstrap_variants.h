#include <random>
#include "multinomial_sampler.h"
#include <unordered_map>

#ifndef IVAR_BOOTSTRAP_VARIANTS_H
#define IVAR_BOOTSTRAP_VARIANTS_H

class bootstrap_variants {
 private:
  unsigned int seed;
  std::mt19937 rng;

  std::vector<uint32_t> sites;
  std::vector<uint32_t> depths;
  std::vector<uint32_t> total_depths;

  std::vector<uint32_t> unique_sites;
  std::vector<uint32_t> unique_sites_total_depths;

  multinomial_sampler* site_sampler;

  std::unordered_map<uint32_t, std::vector<uint32_t>> site_index;

  void build_site_index();

 public:
  bootstrap_variants(const std::vector<uint32_t> &sites_, const std::vector<uint32_t> &depths_, const std::vector<uint32_t> &total_depths_, unsigned int seed = std::random_device{}());
  ~bootstrap_variants() {
    delete site_sampler;
  }

  void sample_sites(std::vector<uint32_t> &sampled_sites, uint32_t n_sites);
  void sample(std::vector<uint32_t> &sampled_variants, std::vector<uint32_t> &sampled_depths, std::vector<double> &sampled_frequencies, uint32_t n_sites);
};

#endif
