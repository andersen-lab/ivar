#include "distribution_sampler.h"
#include <random>
#include <stdexcept>

#ifndef IVAR_MULTINOMIAL_SAMPLER_H
#define IVAR_MULTINOMIAL_SAMPLER_H

class multinomial_sampler: public distribution_sampler  {
 private:
  std::vector<uint32_t> weights;
  unsigned int seed;
  std::mt19937 rng;
  std::discrete_distribution<uint32_t> dist;

 public:
  multinomial_sampler(const std::vector<uint32_t> weights_, unsigned int seed = std::random_device{}()): weights(weights_), seed(seed), rng(seed), dist(weights.begin(), weights.end()) {};

  int sample() override;
  void sample(std::vector<uint32_t> &out, uint32_t n) override;
  void set_seed(unsigned int seed_) override;
  void sample(std::vector<double> &out, uint32_t n) override {
    throw std::runtime_error("multinomial_sampler::sample: not implemented for double");
  }
};

#endif  // IVAR_MULTINOMIAL_SAMPLER_H
