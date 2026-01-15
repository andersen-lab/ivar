#include "../src/gaussian_sampler.h"
#include <iostream>

int main() {
  std::vector<double> samples;
  gaussian_sampler sampler(0.0, 1.0, 42);

  sampler.sample(samples, 5);

  if (samples.size() != 5)
    return -1;

  sampler.sample(samples, 1000);

  if (samples.size() != 1000)
    return -1;

  std::cerr << "Samples: ";
  for (double s : samples)
    std::cerr << s << ", ";

  return 0;
}