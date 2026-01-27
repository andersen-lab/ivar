#include "multinomial_sampler.h"

int multinomial_sampler::sample() {
  return this->dist(this->rng);
}

void multinomial_sampler::sample(std::vector<uint32_t> &out, uint32_t n) {
  out.clear();
  out.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    out.push_back(this->dist(this->rng));
  }
}

void multinomial_sampler::set_seed(uint32_t seed_) {
  this->rng.seed(seed_);
}