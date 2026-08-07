#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

int main() {
  int num_tests = 0;
  int success = 0;

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
