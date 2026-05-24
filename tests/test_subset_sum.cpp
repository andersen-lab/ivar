#include <iostream>
#include <vector>

#include "../src/solve_clustering.h"

int main() {
  int retval = 0;
  double error = 0.05;
  std::vector<double> means = {0.898281, 0.421816, 0.541234, 0.144119, 0.297895}; // 0.596343

  std::vector<std::vector<double>> results;
  subset_sum(means, results, error);

  return 0;
}
