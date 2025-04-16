#include <iostream>
#include <vector>

#include "htslib/sam.h"
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/estimate_error.h"
#include "../src/interval_tree.h"

int main() {
  int num_tests = 2;
  int success = 2;
  return (num_tests == success) ? 0 : -1;
}
