#include "solve_clustering.h"
#include "call_consensus_clustering.h"
#include "genomic_position.h"
#include "gmm.h"
#include "saga.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <limits>
#include <set>

const double PEAK_COLLISION_TOLERANCE = 0.03;

double find_nearest_distance(const std::vector<double> all_sums, double value) {
    double min_distance = std::numeric_limits<double>::max();
    for (auto num : all_sums) {
        double distance = std::abs(num - value);
        if (distance < min_distance) {
            min_distance = distance;
        }
    }
    return min_distance;
}

bool account_for_inferred_means(std::vector<double> possible_solution, std::vector<double> means, double total, double error){
  bool valid = true;
  std::vector<double> current;
  std::vector<std::vector<double>> results;

  find_combinations(possible_solution, 0, current, results, 0);

  std::vector<double> all_sums;
  for(auto result : results){ 
    if(result.size() == 1) continue;
    double sum = std::accumulate(result.begin(), result.end(), 0.0);
    all_sums.push_back(sum);
  }

  //check if all means can be accounted for
  for(auto mean : means){
    bool found = std::find(possible_solution.begin(), possible_solution.end(), mean) != possible_solution.end();
    if(found) continue;
    double dist = find_nearest_distance(all_sums, mean);
    if(dist > error){
      valid = false;
      break;
    }
  }
  return(valid);
}

bool within_error_range(std::vector<double> values, double target, double error){
  //test if the sum of the vector equals the target value within some error
  double sum = std::accumulate(values.begin(), values.end(), 0.0);
  if(sum < target+(error) && sum > target-(error)){
    return(true);
  } else{
    return(false);
  }
}

std::vector<std::vector<double>> find_subsets(std::vector<double> means, double target, double error){
  //first we find all the possible combinations
  std::vector<double> current;
  std::vector<std::vector<double>> results;
  find_combinations(means, 0, current, results, 0);
  std::vector<std::vector<double>> valid_combinations;

  for(uint32_t i=0; i < results.size(); i++){
    bool in_range = within_error_range(results[i], target, error);
    if(in_range){
      valid_combinations.push_back(results[i]);
    }
  }
  return(valid_combinations);
}

void find_combinations(std::vector<double> means, uint32_t index, std::vector<double> &current, std::vector<std::vector<double>> &results, double error){
  if (!current.empty()){
    results.push_back(current);
  }
  for (uint32_t i = index; i < means.size(); ++i) {
    if(means[i] < error) continue;
    current.push_back(means[i]);
    find_combinations(means, i+1, current, results, error);
    current.pop_back();
  }
}

std::vector<std::vector<double>> find_solutions(std::vector<double> means, double error){
  std::vector<double> current;
  std::vector<std::vector<double>> results;

  find_combinations(means, 0, current, results, 0);
  std::sort(results.begin(), results.end());
  results.erase(std::unique(results.begin(), results.end()), results.end());

  std::vector<std::vector<double>> final_results;

  //constrain that the solutions must add to 1
  for(uint32_t i=0; i < results.size(); i++){
    bool keep = within_error_range(results[i], 1, error);
    if(keep){
      final_results.push_back(results[i]);
    }
  }
  return(final_results);
}

std::vector<std::vector<uint32_t>> find_combination_peaks(std::vector<double> solution, std::vector<double> means, std::vector<double> &unresolved, double error, std::vector<std::vector<uint32_t>> &ambiguous_indexes){

  std::vector<std::vector<uint32_t>> cluster_indexes(means.size());
  ambiguous_indexes.assign(means.size(), {});
  std::vector<double> current;
  std::vector<std::vector<double>> results;
  std::vector<double> totals;

  find_combinations(solution, 0, current, results, error);
  for(uint32_t i=0; i < results.size(); i++){
    double sum = std::accumulate(results[i].begin(), results[i].end(), 0.0);
    totals.push_back(sum);
  }

  //given a solution and the means, map each peak to every combination of solution
  //values ("explanation") whose sum lands within PEAK_COLLISION_TOLERANCE of it
  for(uint32_t i=0; i < means.size(); i++){
    double target = means[i];
    auto it = std::find(solution.begin(), solution.end(), target);

    std::vector<std::vector<uint32_t>> explanations;
    if(it != solution.end()){
      explanations.push_back({i}); //the mean is itself a standalone population
    }
    for(uint32_t r=0; r < results.size(); r++){
      if(results[r].size() < 2) continue; //singletons are covered above
      if(std::abs(totals[r] - target) >= PEAK_COLLISION_TOLERANCE) continue;
      std::vector<uint32_t> members;
      for(auto x : results[r]){
        auto it2 = std::find(std::begin(means), std::end(means), x);
        members.push_back((uint32_t)std::distance(std::begin(means), it2));
      }
      explanations.push_back(members);
    }

    if(explanations.empty()){
      //no exact/combination explanation - fall back to nearest total
      auto nearest = std::min_element(totals.begin(), totals.end(), [target](double a, double b) {return std::abs(a - target) < std::abs(b - target);});
      uint32_t index = std::distance(totals.begin(), nearest);
      std::vector<uint32_t> members;
      for(auto x : results[index]){
        auto it2 = std::find(std::begin(means), std::end(means), x);
        members.push_back((uint32_t)std::distance(std::begin(means), it2));
      }
      explanations.push_back(members);
    }

    if(explanations.size() == 1){
      cluster_indexes[i] = explanations[0];
      continue;
    }

    //multiple explanations: confirmed = present in every explanation,
    //ambiguous = present in some but not all
    std::set<uint32_t> confirmed(explanations[0].begin(), explanations[0].end());
    std::set<uint32_t> unioned(explanations[0].begin(), explanations[0].end());
    for(uint32_t e=1; e < explanations.size(); e++){
      std::set<uint32_t> current_set(explanations[e].begin(), explanations[e].end());
      std::set<uint32_t> new_confirmed;
      std::set_intersection(confirmed.begin(), confirmed.end(), current_set.begin(), current_set.end(), std::inserter(new_confirmed, new_confirmed.begin()));
      confirmed = new_confirmed;
      unioned.insert(current_set.begin(), current_set.end());
    }
    cluster_indexes[i] = std::vector<uint32_t>(confirmed.begin(), confirmed.end());
    std::set<uint32_t> ambig;
    std::set_difference(unioned.begin(), unioned.end(), confirmed.begin(), confirmed.end(), std::inserter(ambig, ambig.begin()));
    ambiguous_indexes[i] = std::vector<uint32_t>(ambig.begin(), ambig.end());
    unresolved.push_back(target);
  }
  return(cluster_indexes);
}

bool is_boundary_rescue(const std::vector<double>& means) {
  if (means.size() != 2) return false;
  const double BOUNDARY_TOL = 0.01;
  int boundary_idx = -1, major_idx = -1;
  for (int i = 0; i < 2; ++i) {
    if (std::abs(means[i] - 0.03) < BOUNDARY_TOL)
      boundary_idx = i;
    else if (means[i] > 0.70 && std::abs(means[i] - 0.97) > BOUNDARY_TOL)
      major_idx = i;
  }
  return (boundary_idx >= 0 && major_idx >= 0);
}

bool subset_sum(std::vector<double> &means, std::vector<std::vector<double>> &solution_sets, const double unit_sum_error){
  double combination_error = 0.05;

  //gives all solutions that sum to 1
  std::vector<std::vector<double>> solutions = find_solutions(means, unit_sum_error);
  if(solutions.size() == 0){
    std::cerr << "no solutions found" << std::endl;
    // Rescue: 2 effective means, one absorbed at the 0.03 lower boundary and
    // one intermediate (>0.50). Infer the minor frequency as the complement
    // of the major mean and treat as a valid 2-population solution.
    if (is_boundary_rescue(means)) {
      const double BOUNDARY_TOL = 0.01;
      int boundary_idx = -1, major_idx = -1;
      for (int i = 0; i < 2; ++i) {
        if (std::abs(means[i] - 0.03) < BOUNDARY_TOL) boundary_idx = i;
        else if (means[i] > 0.70) major_idx = i;
      }
      double inferred_minor = 1.0 - means[major_idx];
      std::cerr << "Rescue: lower boundary absorbed minor population; "
                << "inferred complement = " << inferred_minor << "\n";
      means[boundary_idx] = inferred_minor;
      solution_sets.push_back({means[major_idx], inferred_minor});
      return true;
    }
    return(false);
  }

  std::vector<double> non_subset_means;
  for(uint32_t i=0; i < means.size(); i++){
    std::vector<std::vector<double>> tmp = find_subsets(means, means[i], unit_sum_error);
    if(tmp.size() <= 1){
      non_subset_means.push_back(means[i]);
    }
  }

  //reduce solution space to things that contain the non subset peaks
  std::vector<std::vector<double>> realistic_solutions;
  for(uint32_t i=0; i < solutions.size(); i++){
      std::vector<double> tmp = solutions[i];
      bool found = std::all_of(non_subset_means.begin(), non_subset_means.end(), [&tmp](double value) {return std::find(tmp.begin(), tmp.end(), value) != tmp.end();});
      if(found){
        realistic_solutions.push_back(solutions[i]);
      }
  }

  if(realistic_solutions.size() == 0){
    std::cerr << "no realistic solutions found" << std::endl;
    return(false);
  }

  for(uint32_t i=0; i < realistic_solutions.size(); i++){
    bool keep = account_for_inferred_means(realistic_solutions[i], means, 1, combination_error);
    if(keep){
      solution_sets.push_back(realistic_solutions[i]);
    }
  }
  if(solution_sets.size() == 0){
    return(false);
  } else {
    return(true);
  }

}

void overwrite_cluster_assigned(std::vector<variant> &variants, 
                                std::vector<double> eff_means, 
                                std::vector<double> means){

  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].half_normal_upper || variants[i].half_normal_lower) continue;
    uint32_t cluster_assigned = variants[i].cluster_assigned;
    double mean = means[cluster_assigned];
    // Use nearest-match instead of exact equality to handle floating-point
    // precision differences between model_means and eff_means.
    uint32_t best_index = 0;
    double best_dist = std::abs(eff_means[0] - mean);
    for(uint32_t j = 1; j < eff_means.size(); j++){
      double dist = std::abs(eff_means[j] - mean);
      if(dist < best_dist){
        best_dist = dist;
        best_index = j;
      }
    }
    variants[i].cluster_assigned = best_index;
  }                                 

}

void assign_variants_solution(std::vector<double> solution,
                              std::vector<variant> &variants,
                              std::vector<double> means,
                              double threshold){
                    
  double error = 0.05;
  std::vector<double> unresolved;
  std::vector<std::vector<uint32_t>> ambiguous_groups;
  std::vector<std::vector<uint32_t>> cluster_groups = find_combination_peaks(solution, means, unresolved, error, ambiguous_groups);
  std::vector<std::vector<uint32_t>> inverse_groups(means.size());
  std::vector<std::vector<uint32_t>> ambiguous_inverse_groups(means.size());
  for(uint32_t i=0; i < cluster_groups.size(); i++){
    for(uint32_t j=0; j < cluster_groups[i].size(); j++){
      inverse_groups[cluster_groups[i][j]].push_back(i);
    }
    for(uint32_t j=0; j < ambiguous_groups[i].size(); j++){
      ambiguous_inverse_groups[ambiguous_groups[i][j]].push_back(i);
    }
  }

  //check if the variant corresponds to an unresolved cluster
  for(uint32_t i=0; i < variants.size(); i++){
    auto it = std::find(unresolved.begin(), unresolved.end(), means[variants[i].cluster_assigned]);
    if(it != unresolved.end()){
      variants[i].resolved = false;
    }
  }
  //for 100% cases assign all consensus genomes to the variant
  std::vector<uint32_t> all_genomes;
  for(uint32_t i=0; i < solution.size(); i++){
    all_genomes.push_back(i);
  }

  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].half_normal_upper){
      variants[i].consensus_numbers = all_genomes;
    }
  }

  //assign the number of the consensus genome
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].half_normal_upper || variants[i].half_normal_lower) continue;

    //assign to every cluster that's still plausible given the posterior
    //probabilities, not just the single argmax cluster_assigned - a variant
    //within threshold of the best probability is genuinely ambiguous between
    //those clusters, so it should reach every genome any of them map to.
    std::vector<uint32_t> top_clusters;
    if(!variants[i].probabilities.empty()){
      double best_probability = *std::max_element(variants[i].probabilities.begin(), variants[i].probabilities.end());
      for(uint32_t c=0; c < variants[i].probabilities.size(); c++){
        if(variants[i].probabilities[c] > best_probability / threshold){
          top_clusters.push_back(c);
        }
      }
    } else {
      top_clusters.push_back(variants[i].cluster_assigned);
    }

    for(uint32_t j=0; j < inverse_groups.size(); j++){
      //check to make sure you're lookin at a group that's part of the solution
      auto mit = std::find(solution.begin(), solution.end(), means[j]);
      if(mit == solution.end()) continue;

      //assign the point to all applicable groups if any top cluster is a member
      bool matches = false;
      for(auto c : top_clusters){
        if(std::find(inverse_groups[j].begin(), inverse_groups[j].end(), c) != inverse_groups[j].end()){
          matches = true;
          break;
        }
      }
      if(matches){
        variants[i].consensus_numbers.push_back((uint32_t)std::distance(solution.begin(), mit));
      }
    }

    for(uint32_t j=0; j < ambiguous_inverse_groups.size(); j++){
      auto mit = std::find(solution.begin(), solution.end(), means[j]);
      if(mit == solution.end()) continue;

      bool matches = false;
      for(auto c : top_clusters){
        if(std::find(ambiguous_inverse_groups[j].begin(), ambiguous_inverse_groups[j].end(), c) != ambiguous_inverse_groups[j].end()){
          matches = true;
          break;
        }
      }
      if(matches){
        variants[i].ambiguous_numbers.push_back((uint32_t)std::distance(solution.begin(), mit));
      }
    }
  }
}
