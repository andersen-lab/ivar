#include <vector>
#include <cstdint>

#ifndef IVAR_DISTRIBUTION_SAMPLER_H
#define IVAR_DISTRIBUTION_SAMPLER_H

class distribution_sampler {
 public:
  virtual ~distribution_sampler() = default;

  virtual int sample() = 0;
  virtual void sample(std::vector<double> &out, uint32_t n) = 0;
  virtual void sample(std::vector<uint32_t> &out, uint32_t n) = 0;
  virtual void set_seed(uint32_t seed) = 0;
};

#endif  // IVAR_DISTRIBUTION_SAMPLER_H
