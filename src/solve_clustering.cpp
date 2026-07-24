#include "solve_clustering.h"
#include "call_consensus_clustering.h"
#include "genomic_position.h"
#include "gmm.h"
#include "saga.h"
#include <ostream>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <limits>

std::vector<uint32_t> find_missing_indexes(const std::vector<uint32_t>& tmp, const std::vector<uint32_t>& amplicons_to_mask) {
  std::unordered_set<uint32_t> mask_set(amplicons_to_mask.begin(), amplicons_to_mask.end());
  std::vector<uint32_t> missing_indexes;
  for (uint32_t i = 0; i < tmp.size(); ++i) {
    if (mask_set.find(tmp[i]) == mask_set.end()) {
      missing_indexes.push_back(i);
    }
  }
  return missing_indexes;
}

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

bool account_peaks(std::vector<double> possible_solution, std::vector<double> means, double total, double error){
  bool valid = true;
  std::vector<double> current;
  std::vector<std::vector<double>> results;

  find_combinations(possible_solution, 0, current, results, 0);

  std::vector<double> all_sums;
  for(auto result : results){ 
    if(result.size() == 1) continue;
    double sum = std::accumulate(result.begin(), result.end(), 0.0f);
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
  double sum = std::accumulate(values.begin(), values.end(), 0.0f);
  if(sum < target+(error) && sum > target-(error)){
    return(true);
  } else{
    return(false);
  }
}

std::vector<std::vector<double>> find_subsets_with_error(std::vector<double> means, double target, double error){
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

std::vector<std::vector<double>> frequency_pair_finder(std::vector<variant> variants, std::vector<double> means){
  std::vector<std::vector<double>> pairs;
  std::vector<uint32_t> track_positions;

  for(uint32_t i=0; i < variants.size(); i++){
    if(!variants[i].depth_flag && !variants[i].qual_flag && !variants[i].outside_freq_range && variants[i].cluster_assigned != -1){
      auto it = std::find(track_positions.begin(), track_positions.end(), variants[i].position);
      //found
      if(it != track_positions.end()){
        size_t index = std::distance(track_positions.begin(), it);
        pairs[index].push_back(means[variants[i].cluster_assigned]);
      } else{
        std::vector<double> tmp = {means[variants[i].cluster_assigned]};
        pairs.push_back(tmp);
        track_positions.push_back(variants[i].position);
      }
    }
  }

  return(pairs);
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

std::vector<std::vector<uint32_t>> find_combination_peaks(std::vector<double> solution, std::vector<double> means, std::vector<double> &unresolved, double error){

  std::vector<std::vector<uint32_t>> cluster_indexes(means.size());
  std::vector<double> current;
  std::vector<std::vector<double>> results;
  std::vector<double> totals;

  find_combinations(solution, 0, current, results, error);
  for(uint32_t i=0; i < results.size(); i++){
    double sum = std::accumulate(results[i].begin(), results[i].end(), 0.0f);
    totals.push_back(sum);
  }

  //given a solution and the means, map each cluster to the cluster it contains
  for(uint32_t i=0; i < means.size(); i++){
    double target = means[i];
    auto it = std::find(solution.begin(), solution.end(), target);

    //the mean is part of the solution
    if(it != solution.end()){
        cluster_indexes[i].push_back(i);
        std::vector<double> distances(totals.size());
        std::transform(totals.begin(), totals.end(), distances.begin(), [target](double num) { return std::abs(target - num); });
        uint32_t count = 0;
        //this checks the distances from the mean to all other possible peaks
        for(uint32_t d=0; d < distances.size(); d++){
          if(distances[d] < 0.03 && distances[d] != 0){
            count += 1;
          }
        }
        if(count > 1) unresolved.push_back(target);

    } else {
      //the problem with this is that it looks at the min but not if two overlapping peaks occur
      auto it = std::min_element(totals.begin(), totals.end(), [target](double a, double b) {return std::abs(a - target) < std::abs(b - target);});

      std::vector<double> distances(totals.size());
      std::transform(totals.begin(), totals.end(), distances.begin(), [target](double num) { return std::abs(target - num); });
      uint32_t count = 0;
      for(uint32_t d=0; d < distances.size(); d++){
        if(distances[d] < 0.03) count += 1;
      }
      uint32_t index = std::distance(totals.begin(), it);
      for(auto x : results[index]){
        auto it2 = std::find(std::begin(means), std::end(means), x);
        uint32_t index2 = std::distance(std::begin(means), it2);
        cluster_indexes[i].push_back(index2);
      }
      if(count > 1) unresolved.push_back(means[i]);
    }
  }
  /*for(uint32_t i=0; i < cluster_indexes.size(); i++){
    for(uint32_t j=0; j < cluster_indexes[i].size(); j++){
      std::cerr << cluster_indexes[i][j] << " ";
    }
    std::cerr << "\n";
  }*/
  //for(auto u : unresolved) std::cerr << u << std::endl;
  return(cluster_indexes);
}

std::vector<std::vector<double>> deduplicate_solutions(std::vector<std::vector<double>> vectors){
  std::vector<std::vector<double>> solutions;
  for(uint32_t i=0; i < vectors.size(); i++){
    if(i == 0){
      solutions.push_back(vectors[i]);
      continue;
    }
    bool add = true;
    for(uint32_t j=0; j < solutions.size(); j++){
      bool same = std::equal(solutions[j].begin(), solutions[j].end(), vectors[i].begin());
      if(same && (solutions[j].size() == vectors[i].size())){
        add = false;
        break;
      }
    }
    if(add){
      solutions.push_back(vectors[i]);
    }
  }
  return(solutions);
}

bool is_boundary_rescue(const std::vector<double>& means) {
  if (means.size() != 2) return false;
  const double BOUNDARY_TOL = 0.01;
  int boundary_idx = -1, major_idx = -1;
  for (int i = 0; i < 2; ++i) {
    if (std::abs(means[i] - 0.03) < BOUNDARY_TOL)
      boundary_idx = i;
    else if (means[i] > 0.50 && std::abs(means[i] - 0.97) > BOUNDARY_TOL)
      major_idx = i;
  }
  return (boundary_idx >= 0 && major_idx >= 0);
}

bool subset_sum(std::vector<double> &means, std::vector<std::vector<double>> &solution_sets, const double error){
  double combination_error = 0.05;
  std::cerr << "in subset sum" << std::endl;
  //gives all solutions that sum to 1
  std::vector<std::vector<double>> solutions = find_solutions(means, error);
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
        else if (means[i] > 0.50) major_idx = i;
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
    std::vector<std::vector<double>> tmp = find_subsets_with_error(means, means[i], error);
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
    bool keep = account_peaks(realistic_solutions[i], means, 1, combination_error);
    std::cerr << "keep " << keep << std::endl;
    for(auto s : realistic_solutions[i]){
      std::cerr << s << " ";
    }
    std::cerr << "\n";
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
    if(variants[i].position == 11000){
      std::cerr << "position " << variants[i].position << " freq " << variants[i].freq << " gapped freq " << variants[i].gapped_freq << std::endl;
      std::cerr << "base " << mean << " cluster assigned " << cluster_assigned << std::endl;
      std::cerr << "eff mean index " << best_index << " effective mean " << eff_means[best_index] << std::endl;
      for(auto h : means){
        std::cerr << h << " ";
      }
      std::cerr << "\n";
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
  std::vector<std::vector<uint32_t>> cluster_groups = find_combination_peaks(solution, means, unresolved, error);
  std::vector<std::vector<uint32_t>> inverse_groups(means.size());
  for(uint32_t i=0; i < cluster_groups.size(); i++){
    for(uint32_t j=0; j < cluster_groups[i].size(); j++){
      inverse_groups[cluster_groups[i][j]].push_back(i);
    }
  }

  double largest = *std::max_element(solution.begin(), solution.end());
  //define the clusters which contain the majority population
  std::vector<std::vector<double>> possible_clusters;
  std::vector<double> current;
  find_combinations(solution, 0, current, possible_clusters, 0);
  std::vector<double> expected_clusters;
  for(uint32_t i=0; i < possible_clusters.size(); i++){
    bool keep = false;
    for(uint32_t j=0; j < possible_clusters[i].size(); j++){
      if(possible_clusters[i][j] == largest){
        keep = true;
        break;
      }
    }
    if(keep){
      double sum = std::accumulate(possible_clusters[i].begin(), possible_clusters[i].end(), 0.0f);
      expected_clusters.push_back(sum);
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
        variants[i].consensus_numbers.push_back(j);
      }
    }
  }
}
