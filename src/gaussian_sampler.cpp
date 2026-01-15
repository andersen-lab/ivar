#include "gaussian_sampler.h"

int gaussian_sampler::sample() {
    return this->dist(this->rng);
}

void gaussian_sampler::sample(std::vector<double> &out, uint32_t n) {
  out.clear();
  out.reserve(n);

  for (int i = 0; i < n; ++i) {
    out.push_back(this->dist(this->rng));
  }
}

void gaussian_sampler::set_seed(uint32_t seed) {
  this->rng.seed(seed);
}

