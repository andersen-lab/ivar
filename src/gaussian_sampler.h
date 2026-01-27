#include "distribution_sampler.h"
#include <random>
#include <stdexcept>

#ifndef IVAR_GAUSSIAN_SAMPLER_H
#define IVAR_GAUSSIAN_SAMPLER_H

class gaussian_sampler: public distribution_sampler {
 private:
  double mean;
  double stddev;
  std::mt19937 rng;
  std::normal_distribution<double> dist;

 public:
  gaussian_sampler(double mean, double stddev, uint32_t seed = std::random_device{}()): mean(mean), stddev(stddev), rng(seed), dist(mean, stddev) {};

  int sample() override;
  void sample(std::vector<double> &out, uint32_t n) override;
  void sample(std::vector<uint32_t> &out, uint32_t n) override {
    throw std::runtime_error("gaussian_sampler::sample: not implemented for uint32_t");
  }
  void set_seed(uint32_t seed) override;
};

#endif  // IVAR_GAUSSIAN_SAMPLER_H
