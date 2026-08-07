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

  const double error = 0.10;

  //TEST 1 - Since our error is 0.10 we should be able to find a solution for means 0.11 and 0.80
  std::vector<double> means = {0.11, 0.80};
  std::vector<std::vector<double>> solution_sets;
  bool solution_status = subset_sum(means, solution_sets, error);
  num_tests++;
  if(solution_sets.size() == 1) {
    success++;
  } else {
    std::cerr << "TEST 1 failed, found " << solution_sets.size() << " solutions" << std::endl;
  }

  //TEST 2 - Since our error is 0.10 we should NOT be able to find a solution for means 0.10 and 0.80
  means.clear();
  solution_sets.clear();
  means = {0.10, 0.80};
  solution_status = subset_sum(means, solution_sets, error);
  num_tests++;
  if(solution_sets.size() == 0) {
    success++;
  } else {
    std::cerr << "TEST 2 failed, found " << solution_sets.size() << " solutions" << std::endl;
  }

  //TEST 3 - We have a 0.1, 0.2, 0.3, 0.4 mixture with mutations shared by 0.2 and 0.3
  // this should be solveable
  means.clear();
  solution_sets.clear();
  means = {0.10, 0.20, 0.30, 0.40, 0.50};
  solution_status = subset_sum(means, solution_sets, error);
  num_tests++;
  if(solution_sets.size() == 1) {
    success++;
  } else {
    std::cerr << "TEST 3 failed, found " << solution_sets.size() << " solutions" << std::endl;
  }

  //TEST 4 - The boundary rescue should infer a lower population if it gets wrapped into the lower half normal
  means.clear();
  solution_sets.clear();
  means = {0.03, 0.85};
  solution_status = subset_sum(means, solution_sets, error);
  num_tests++;
  if(solution_sets.size() == 1) {
    success++;
  } else {
    std::cerr << "TEST 4 failed, found " << solution_sets.size() << " solutions" << std::endl;
  }
  
  //TEST 5 - This should have too many solutions and not be solveable
  means.clear();
  solution_sets.clear();
  means = {0.10, 0.12, 0.15, 0.18, 0.20, 0.30, 0.50};
  solution_status = subset_sum(means, solution_sets, error);
  num_tests++;
  if(solution_sets.size() == 3) {
    success++;
  } else {
    std::cerr << "TEST 5 failed, found " << solution_sets.size() << " solutions" << std::endl;
  }
  
  //TEST 6 - This should have no solutions and not be solveable
  means.clear();
  solution_sets.clear();
  means = {0.10, 0.30, 0.50, 0.70};
  solution_status = subset_sum(means, solution_sets, error);
  num_tests++;
  if(solution_sets.size() == 0) {
    success++;
  } else {
    std::cerr << "TEST 6 failed, found " << solution_sets.size() << " solutions" << std::endl;
  }
  
  /*std::cerr << solution_sets.size() << " solutions found for means: ";
  for(auto solution : solution_sets){
    for(auto s : solution){
      std::cerr << s << " ";
    }
    std::cerr << "\n";
  }*/
  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
