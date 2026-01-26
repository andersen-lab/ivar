#include "distribution_sampler.h"
#include <random>

#ifndef IVAR_MULTINOMIAL_SAMPLER_H
#define IVAR_MULTINOMIAL_SAMPLER_H

class multinomial_sampler: public distribution_sampler  {
 private:
  std::vector<double> weights;
  unsigned int seed;
  std::mt19937 rng;
  std::discrete_distribution<int> dist;

 public:
  multinomial_sampler(const std::vector<double> weights_, uint32_t seed = std::random_device{}()): weights(weights_), seed(seed), rng(seed), dist(weights.begin(), weights.end()) {};

  int sample() override;
  void sample(std::vector<double> &out, uint32_t n) override;
  void set_seed(uint32_t seed_) override;
};

#endif  // IVAR_MULTINOMIAL_SAMPLER_H
