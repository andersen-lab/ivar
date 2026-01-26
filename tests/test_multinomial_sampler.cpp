#include "../src/multinomial_sampler.h"
#include <iostream>

int main() {
  std::vector<double> samples;
  std::vector<double> weights = {0.1,0.6,0.3};
  multinomial_sampler sampler(weights, 112358);

  int idx = sampler.sample();
  if(idx < 0 || idx >= weights.size())
    return -1;

  sampler.sample(samples, 100);
  if (samples.size() != 100)
    return -1;

  return 0;
}