#include "solve_clustering.h"
#include "gmm.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <limits>
#include <set>

constexpr double subset_sum_solver::SUBSET_SUM_ERROR;
constexpr double subset_sum_solver::MAJOR_MEAN_FLOOR;
constexpr double subset_sum_solver::UNIT_SUM_ERROR;
constexpr double variant_assigner::PEAK_COLLISION_TOLERANCE;
constexpr double variant_assigner::COMBINATION_MIN_MEAN;

//Shared by both classes, so it lives at file scope with internal linkage rather
//than on either one. Every non-empty subset of `values`, skipping entries below
//`min_value` - note that is an element floor, not a tolerance.
static void enumerate_recursive(const std::vector<double> &values, uint32_t index, std::vector<double> &current, std::vector<std::vector<double>> &results, double min_value){
  if (!current.empty()){
    results.push_back(current);
  }
  for (uint32_t i = index; i < values.size(); ++i) {
    if(values[i] < min_value) continue;
    current.push_back(values[i]);
    enumerate_recursive(values, i+1, current, results, min_value);
    current.pop_back();
  }
}

static std::vector<std::vector<double>> enumerate_combinations(const std::vector<double> &values, double min_value){
  std::vector<double> current;
  std::vector<std::vector<double>> results;
  enumerate_recursive(values, 0, current, results, min_value);
  return results;
}

// ---------------------------------------------------------------- subset_sum_solver

subset_sum_solver::subset_sum_solver(std::vector<double> means, double unit_sum_error, double invariant_threshold)
  : means(means),
    unit_sum_error(unit_sum_error),
    invariant_mean(invariant_threshold),
    boundary_mean(1.0 - invariant_threshold) {}

double subset_sum_solver::find_nearest_distance(const std::vector<double> &all_sums, double value) {
  double min_distance = std::numeric_limits<double>::max();
  for (auto num : all_sums) {
    double distance = std::abs(num - value);
    if (distance < min_distance) {
      min_distance = distance;
    }
  }
  return min_distance;
}

bool subset_sum_solver::within_error_range(const std::vector<double> &values, double target, double error){
  //test if the sum of the vector equals the target value within some error
  //not equivalent to abs(sum-target) < error in floating point, see TEST 2
  double sum = std::accumulate(values.begin(), values.end(), 0.0);
  return sum < target+(error) && sum > target-(error);
}

bool subset_sum_solver::account_for_inferred_means(const std::vector<double> &candidate) const {
  std::vector<std::vector<double>> results = enumerate_combinations(candidate, 0);

  std::vector<double> all_sums;
  for(auto result : results){
    if(result.size() == 1) continue;
    double sum = std::accumulate(result.begin(), result.end(), 0.0);
    all_sums.push_back(sum);
  }

  //check if all means can be accounted for
  for(auto mean : means){
    bool found = std::find(candidate.begin(), candidate.end(), mean) != candidate.end();
    if(found) continue;
    double dist = find_nearest_distance(all_sums, mean);
    //exclusive, matching within_error_range
    if(dist >= SUBSET_SUM_ERROR){
      return false;
    }
  }
  return true;
}

std::vector<std::vector<double>> subset_sum_solver::find_subsets(double target) const {
  std::vector<std::vector<double>> results = enumerate_combinations(means, 0);
  std::vector<std::vector<double>> valid_combinations;

  for(uint32_t i=0; i < results.size(); i++){
    if(within_error_range(results[i], target, SUBSET_SUM_ERROR)){
      valid_combinations.push_back(results[i]);
    }
  }
  return valid_combinations;
}

std::vector<std::vector<double>> subset_sum_solver::find_solutions() const {
  std::vector<std::vector<double>> results = enumerate_combinations(means, 0);
  //duplicate values in `means` produce equal subsets, so this is not a no-op
  std::sort(results.begin(), results.end());
  results.erase(std::unique(results.begin(), results.end()), results.end());

  std::vector<std::vector<double>> final_results;

  //constrain that the solutions must add to 1
  for(uint32_t i=0; i < results.size(); i++){
    if(within_error_range(results[i], 1, unit_sum_error)){
      final_results.push_back(results[i]);
    }
  }
  return final_results;
}

//returns the indices it found rather than discarding them, so detection and
//application cannot drift apart. gmm_1d pins the half normal means to exactly
//1-invariant_threshold and invariant_threshold (gmm_1d.cpp:306) and
//get_effective_means copies them untransformed, so these comparisons are exact
bool subset_sum_solver::find_boundary_rescue(int &boundary_idx, int &major_idx) const {
  boundary_idx = -1;
  major_idx = -1;
  if (means.size() != 2) return false;
  for (int i = 0; i < 2; ++i) {
    if (means[i] == boundary_mean)
      boundary_idx = i;
    else if (means[i] > MAJOR_MEAN_FLOOR && means[i] != invariant_mean)
      major_idx = i;
  }
  return (boundary_idx >= 0 && major_idx >= 0);
}

// Rescue: 2 effective means, one absorbed at the lower half normal boundary
// (1-invariant_threshold) and one major peak above MAJOR_MEAN_FLOOR that is not the
// invariant peak. Infer the minor frequency as the complement of the major mean and
// treat as a valid 2-population solution.
bool subset_sum_solver::apply_boundary_rescue(){
  int boundary_idx = -1, major_idx = -1;
  if(!find_boundary_rescue(boundary_idx, major_idx)) return false;

  double inferred_minor = 1.0 - means[major_idx];
  std::cerr << "Rescue: lower boundary absorbed minor population; "
            << "inferred complement = " << inferred_minor << "\n";
  //the same double lands in both, which downstream exact-equality matching needs
  means[boundary_idx] = inferred_minor;
  solution_sets.push_back({means[major_idx], inferred_minor});
  return true;
}

bool subset_sum_solver::solve(){
  //gives all solutions that sum to 1
  std::vector<std::vector<double>> solutions = find_solutions();
  if(solutions.size() == 0){
    std::cerr << "no solutions found" << std::endl;
    return apply_boundary_rescue();
  }

  std::vector<double> non_subset_means;
  for(uint32_t i=0; i < means.size(); i++){
    if(find_subsets(means[i]).size() <= 1){
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
    return false;
  }

  for(uint32_t i=0; i < realistic_solutions.size(); i++){
    if(account_for_inferred_means(realistic_solutions[i])){
      solution_sets.push_back(realistic_solutions[i]);
    }
  }
  return !solution_sets.empty();
}

// ----------------------------------------------------------------- variant_assigner

variant_assigner::variant_assigner(std::vector<double> solution, std::vector<double> means, double threshold)
  : solution(solution), means(means), threshold(threshold) {}

std::vector<uint32_t> variant_assigner::members_of(const std::vector<double> &subset) const {
  std::vector<uint32_t> members;
  for(auto x : subset){
    auto it = std::find(means.begin(), means.end(), x);
    members.push_back((uint32_t)std::distance(means.begin(), it));
  }
  return members;
}

std::vector<std::vector<uint32_t>> variant_assigner::find_combination_peaks(std::vector<double> &unresolved, std::vector<std::vector<uint32_t>> &ambiguous_indexes) const {
  std::vector<std::vector<uint32_t>> cluster_indexes(means.size());
  ambiguous_indexes.assign(means.size(), {});

  std::vector<std::vector<double>> results = enumerate_combinations(solution, COMBINATION_MIN_MEAN);
  std::vector<double> totals;
  for(uint32_t i=0; i < results.size(); i++){
    totals.push_back(std::accumulate(results[i].begin(), results[i].end(), 0.0));
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
      explanations.push_back(members_of(results[r]));
    }

    if(explanations.empty()){
      //no exact/combination explanation - fall back to nearest total
      auto nearest = std::min_element(totals.begin(), totals.end(), [target](double a, double b) {return std::abs(a - target) < std::abs(b - target);});
      uint32_t index = std::distance(totals.begin(), nearest);
      explanations.push_back(members_of(results[index]));
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
  return cluster_indexes;
}

void variant_assigner::overwrite_cluster_assigned(std::vector<variant> &variants,
                                                  const std::vector<double> &eff_means,
                                                  const std::vector<double> &means){
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

//assign to every cluster that's still plausible given the posterior probabilities,
//not just the single argmax cluster_assigned - a variant within threshold of the best
//probability is genuinely ambiguous between those clusters, so it should reach every
//genome any of them map to.
void variant_assigner::collect_matches(const std::vector<std::vector<uint32_t>> &inverse_groups,
                                       const std::vector<uint32_t> &top_clusters,
                                       std::vector<uint32_t> &out) const {
  for(uint32_t j=0; j < inverse_groups.size(); j++){
    //check to make sure you're lookin at a group that's part of the solution
    auto mit = std::find(solution.begin(), solution.end(), means[j]);
    if(mit == solution.end()) continue;

    bool matches = false;
    for(auto c : top_clusters){
      if(std::find(inverse_groups[j].begin(), inverse_groups[j].end(), c) != inverse_groups[j].end()){
        matches = true;
        break;
      }
    }
    if(matches){
      out.push_back((uint32_t)std::distance(solution.begin(), mit));
    }
  }
}

void variant_assigner::assign(std::vector<variant> &variants) const {
  std::vector<double> unresolved;
  std::vector<std::vector<uint32_t>> ambiguous_groups;
  std::vector<std::vector<uint32_t>> cluster_groups = find_combination_peaks(unresolved, ambiguous_groups);
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

    collect_matches(inverse_groups, top_clusters, variants[i].consensus_numbers);
    collect_matches(ambiguous_inverse_groups, top_clusters, variants[i].ambiguous_numbers);
  }
}
