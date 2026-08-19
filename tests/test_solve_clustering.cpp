#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
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
  
  //TEST 7 - Same 0.1/0.2/0.3/0.4 mixture as TEST 3, but this time exercising
  //assign_variants_solution directly. A variant sitting on the 0.40 peak is
  //explainable two ways - {0.40} alone, or {0.10 + 0.30} - so it should be
  //*ambiguous* across 0.10/0.30/0.40 (none of them confirmed) rather than
  //confidently assigned to the 0.40 genome. A variant on the 0.20 peak has no
  //colliding combination, so it should be confidently assigned as before, with
  //no ambiguity at all.
  {
    std::vector<double> solution = {0.10, 0.20, 0.30, 0.40};
    std::vector<double> peak_means = {0.10, 0.20, 0.30, 0.40};

    variant v_collision;
    v_collision.cluster_assigned = 3; //index of the 0.40 peak in peak_means

    variant v_private;
    v_private.cluster_assigned = 1; //index of the 0.20 peak in peak_means

    std::vector<variant> variants = {v_collision, v_private};
    assign_variants_solution(solution, variants, peak_means, 2.0);

    std::set<uint32_t> collision_confirmed(variants[0].consensus_numbers.begin(), variants[0].consensus_numbers.end());
    std::set<uint32_t> collision_ambiguous(variants[0].ambiguous_numbers.begin(), variants[0].ambiguous_numbers.end());
    std::set<uint32_t> expected_ambiguous = {0, 2, 3}; //0.10, 0.30, 0.40 - not 0.20

    num_tests++;
    if(collision_confirmed.empty() && collision_ambiguous == expected_ambiguous) {
      success++;
    } else {
      std::cerr << "TEST 7a failed, confirmed=[ ";
      for(auto c : collision_confirmed) std::cerr << c << " ";
      std::cerr << "] ambiguous=[ ";
      for(auto a : collision_ambiguous) std::cerr << a << " ";
      std::cerr << "]" << std::endl;
    }

    std::set<uint32_t> private_confirmed(variants[1].consensus_numbers.begin(), variants[1].consensus_numbers.end());
    std::set<uint32_t> private_ambiguous(variants[1].ambiguous_numbers.begin(), variants[1].ambiguous_numbers.end());

    num_tests++;
    if(private_confirmed == std::set<uint32_t>{1} && private_ambiguous.empty()) {
      success++;
    } else {
      std::cerr << "TEST 7b failed, confirmed=[ ";
      for(auto c : private_confirmed) std::cerr << c << " ";
      std::cerr << "] ambiguous=[ ";
      for(auto a : private_ambiguous) std::cerr << a << " ";
      std::cerr << "]" << std::endl;
    }
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
