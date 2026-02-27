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

void call_majority_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth){
  std::cerr << "in majority consensus call" << std::endl;
  uint32_t max_position=0;
  uint32_t min_position = 4294967295U;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
    if(x.position < min_position && x.total_depth > 0){
      min_position = x.position;
    }
  }

  std::vector<std::string> nucs;
  std::vector<double> freqs;
  std::vector<std::string> tmp(max_position, "N");


  for(uint32_t i=1; i <= max_position; i++){
    freqs.clear();
    nucs.clear();
    for(uint32_t j=0; j < variants.size(); j++){
      if(variants[j].position == i){
        if(variants[j].total_depth < min_depth){
          continue;
        }
        nucs.push_back(variants[j].nuc);
        freqs.push_back(variants[j].gapped_freq);
      }
    }
    if(freqs.size() == 0) continue;
    //find the largest frequency
    uint32_t index = std::distance(freqs.begin(), std::max_element(freqs.begin(), freqs.end()));
    if(default_threshold > 0){
      if(freqs[index] >= default_threshold){
        tmp[i-1] = nucs[index];
      }
    } else {
      tmp[i-1] = nucs[index];
    }
  }
  std::cerr << "tmp length " << tmp.size() << std::endl;
  std::string consensus_string = std::accumulate(tmp.begin(), tmp.end(), std::string(""));
  //std::string trimmed_consensus = trim_trailing_ambiguities(consensus_string, max_position);
  std::string next_trimmed_consensus = trim_leading_ambiguities(consensus_string, min_position);
  //write the consensus to file

  std::string consensus_filename = clustering_file + "_threshold.fa";
  std::ofstream file(consensus_filename);
  std::string name = ">"+clustering_file+"_"+std::to_string(default_threshold)+"_threshold";
  file << name << "\n";
  file << next_trimmed_consensus << "\n";
  file.close();
}

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

void modify_variant_masking(std::vector<uint32_t> amplicons_to_mask, std::vector<variant> &variants, std::vector<double> means){
  std::cerr << "modify variant masking" << std::endl;
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].consensus_numbers.size() == 0) continue;
    std::vector<uint32_t> tmp = variants[i].amplicon_numbers;
    std::vector<uint32_t> valid_amplicons = find_missing_indexes(tmp, amplicons_to_mask);
    if(valid_amplicons.size() == 0){
      variants[i].amplicon_masked = true;
    } else if(valid_amplicons.size() == tmp.size()){
      variants[i].amplicon_masked = false;
      variants[i].amplicon_flux = false;
    } else {
      for(auto j : valid_amplicons){
        //std::cerr << variants[i].position << " " << variants[i].cluster_assigned << " " << variants[i].gapped_freq << std::endl;
        variants[i].cluster_assigned = variants[i].consensus_numbers[j];
        variants[i].amplicon_masked = false;
      }
    }
  }
  std::cerr << "end variant masking" << std::endl;
}

bool test_cluster_deviation(double nearest_cluster, double variant_cluster, double std_dev){
  bool fluctuation = false;
  std::vector<double> tmp = {nearest_cluster, variant_cluster};
  double cluster_dev = calculate_standard_deviation(tmp);
  if(std_dev > cluster_dev){
    fluctuation = true;
  }
  return(fluctuation);
}

double find_neighboring_cluster(double freq, uint32_t cluster_assigned, std::vector<double> means){
  //find closest cluster by mean value
  double min_dist = 1;
  uint32_t index = 0;
  for(uint32_t i=0; i < means.size(); i++){
    if(i == cluster_assigned) continue;
    double dist = std::abs(means[i]-freq);
    if(dist < min_dist){
      min_dist = dist;
      index = i;
    }
  }
  return(means[index]);
}

void amplicon_specific_cluster_assignment(std::vector<variant> &variants, gaussian_mixture_model model){
  std::vector<std::vector<double>> prob_matrix;
  std::vector<double> tmp;

  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].freq_numbers.size() < 2) continue;
    if(variants[i].std_dev == 0) continue;
    if(!variants[i].amplicon_flux && !variants[i].amplicon_masked) continue;
    arma::mat final_data = arma::conv_to<arma::rowvec>::from(variants[i].freq_numbers);
    final_data.reshape(1, variants[i].freq_numbers.size());
    tmp.clear();
    prob_matrix.clear();
    for(uint32_t j=0; j < model.n; j++){
      arma::rowvec set_likelihood = model.model.log_p(final_data, j);
      tmp.clear();
      for(uint32_t k=0; k < final_data.n_cols; k++){
        tmp.push_back((double)set_likelihood[k]);
      }
      prob_matrix.push_back(tmp);
    }
    std::vector<std::vector<double>> inverse = transpose_vector(prob_matrix);
    for(uint32_t j=0; j < variants[i].freq_numbers.size(); j++){
      auto max_it = std::max_element(inverse[j].begin(), inverse[j].end());
      uint32_t index = std::distance(inverse[j].begin(), max_it);
      variants[i].freq_assignments.push_back(index);
    }
  }
}

void rewrite_position_masking(std::vector<variant> &variants){
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].freq_numbers.size() < 2) continue;
    if(!variants[i].amplicon_flux) continue;
      bool all_equal = std::all_of(variants[i].freq_assignments.begin(), variants[i].freq_assignments.end(), [&](uint32_t v) {return v == variants[i].freq_assignments[0];});
      if(all_equal) variants[i].amplicon_flux = false;
      else variants[i].amplicon_flux = true;
  }
}

std::vector<uint32_t> rewrite_amplicon_masking(std::vector<variant> variants, std::vector<double> means){
  //stores the numbers of every amplicon where we believe experiences fluctuation that imapcts consensus
  //here we define that as a position where the amplicon is fluctuating and the closest cluster is within a standard deviation
  std::vector<uint32_t> amplicons_to_mask;

  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].amplicon_flux && !variants[i].outside_freq_range){
      //find all clusters not part of the same assignment
      std::vector<double> other_population_clusters;
      other_population_clusters.reserve(means.size());
      for(uint32_t j=0; j< means.size();j++){
        auto it = std::find(variants[i].consensus_numbers.begin(), variants[i].consensus_numbers.end(), j);
        if(it == variants[i].consensus_numbers.end()) other_population_clusters.push_back(means[j]);
      }
      //find the second closest cluster index
      double closest_mean = find_neighboring_cluster(variants[i].gapped_freq, variants[i].cluster_assigned, other_population_clusters);
      //check if the cluster is within the standard dev of the variant
      bool fluctuating = test_cluster_deviation(closest_mean, means[variants[i].cluster_assigned], variants[i].std_dev);
      if(!fluctuating) continue;
      for(auto v : variants[i].amplicon_numbers){
        if(std::find(amplicons_to_mask.begin(), amplicons_to_mask.end(), v) == amplicons_to_mask.end()){
          amplicons_to_mask.push_back(v);
        }
      }
    }
  }
  return(amplicons_to_mask);
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

std::vector<uint32_t> noise_cluster_calculator(gaussian_mixture_model model, double estimated_error){
  std::vector<double> means = model.means;
  std::vector<double> std_devs = model.cluster_std_devs;
  std::vector<uint32_t> noise_indices; 
  for(uint32_t i=0; i < means.size(); i++){
    //if the estimated error is below two standard deviation of the cluster mean
    if(std_devs[i] > 0.05) continue;
    double upper_bound = 1-estimated_error;
    double lower_bound = estimated_error;
    double cluster_lower_edge = means[i] - std_devs[i];
    double cluster_upper_edge = means[i] + std_devs[i];

    if(cluster_lower_edge < lower_bound || cluster_upper_edge > upper_bound){
      //std::cerr << "HERE " << means[i] << " std dev " << std_devs[i] << " estimated error " << estimated_error << " " << cluster_lower_edge << " " << cluster_upper_edge << std::endl;
      noise_indices.push_back(i);
    }
  }
  return(noise_indices);
}

bool subset_sum(std::vector<double> means, std::vector<std::vector<double>> &solution_sets, const double error){
  std::cerr << "in subset sum" << std::endl;
  //gives all solutions that sum to 1
  std::vector<std::vector<double>> solutions = find_solutions(means, error);
  if(solutions.size() == 0){
    std::cerr << "no solutions found" << std::endl;
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
    bool keep = account_peaks(realistic_solutions[i], means, 1, error);
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
    std::cerr << "final no solution sets found" << std::endl;
    return(false);
  } else {
    return(true);
  }

}

std::vector<std::vector<double>> subset_sum(std::vector<double> means){
  double error = 0.05;
  double solution_error = 0.05;

  //determine if any clusters are possible noise
  std::vector<uint32_t> noise_indices;
 
  //filter peaks from means by index
  std::vector<double> filtered_means = means;
  std::vector<double> std_devs;

  std::vector<std::vector<double>> solutions = find_solutions(filtered_means, error);
  //std::cerr << "solution size " << solutions.size() << std::endl;
  //find peaks that can't be a subset of other peaks
  std::vector<double> non_subset_means;
  for(uint32_t i=0; i < filtered_means.size(); i++){
    std::vector<std::vector<double>> tmp = find_subsets_with_error(filtered_means, filtered_means[i], std_devs[i]);
    if(tmp.size() <= 1){
      non_subset_means.push_back(filtered_means[i]);
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
  //std::cerr << "realistic solutions " << realistic_solutions.size() << std::endl;
  //check each solution that every possible peak is accounted for
  std::vector<std::vector<double>> solution_sets;
  for(uint32_t i=0; i < realistic_solutions.size(); i++){
    bool keep = account_peaks(realistic_solutions[i], filtered_means, 1, solution_error);
    if(keep){
      solution_sets.push_back(realistic_solutions[i]);
    }
  }
  //std::cerr << "solution sets " << solution_sets.size() << std::endl;
  return(solution_sets);
}

void solve_clusters(std::vector<variant> &variants, 
                    gaussian_mixture_model model, 
                    double estimated_error, 
                    std::vector<double> &solution, 
                    std::string prefix, 
                    double default_threshold, 
                    uint32_t min_depth){

  double error = 0.05;
  std::vector<double> means = model.means;
  std::cerr << "estimated error " << estimated_error << std::endl;
  std::vector<std::vector<double>> solution_sets;
  for(auto sol : solution_sets){
    for(auto s : sol){
      std::cerr << s << " ";
    }
    std::cerr << "\n";
  }
  bool traditional_majority= false; //if we can't find a solution call a traditional majority consensus
  if(solution_sets.size() == 0){
   std::cerr << "no solution found" << std::endl;
   traditional_majority = true;
  } else if(solution_sets.size() > 1){
    std::cerr << "too many solutions" << std::endl;
    traditional_majority = true;
  } else{
    solution = solution_sets[0];
  }
  if(traditional_majority){
    call_majority_consensus(variants, prefix, default_threshold, min_depth);
    exit(0);
  }
  std::vector<double> unresolved;
  std::vector<std::vector<uint32_t>> cluster_groups = find_combination_peaks(solution, means, unresolved, error);
  std::vector<std::vector<uint32_t>> inverse_groups(means.size());
  for(uint32_t i=0; i < cluster_groups.size(); i++){
    for(uint32_t j=0; j < cluster_groups[i].size(); j++){
      inverse_groups[cluster_groups[i][j]].push_back(i);
    }
  }
  //TESTLINES CHANGE THIS
  std::string solution_string = "[";
  for(uint32_t i=0; i < solution_sets.size(); i++){
    if(i !=0) solution_string += ",";
    solution_string += "[";
      for(uint32_t j=0; j < solution_sets[i].size(); j++){
        if(j != 0){
          solution_string += ",";
        }
        std::string tmp = std::to_string(solution_sets[i][j]);
        solution_string += tmp;
      }
    solution_string +="]";
  }

  solution_string += "]";
  std::string solution_filename = prefix + "_solution.txt";
  std::ofstream file_sol(solution_filename);
  file_sol << "means\n";
  file_sol << solution_string << "\n";
  file_sol.close();
  //exit(0);

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
  //a list of cluster assignments that we assign to consensus
  std::vector<int> major_indexes;
  //index of the "100%" cluster
  for(uint32_t j=0; j < means.size(); j++){
    double tmp = means[j];
    auto closest = *std::min_element(expected_clusters.begin(), expected_clusters.end(), [tmp](double a, double b) {
      return std::abs(a - tmp) < std::abs(b - tmp);
    });
    double diff = std::abs(closest - tmp);
    auto it = std::find(solution.begin(), solution.end(), tmp);
    if((diff < error && it == solution.end()) || tmp == largest){
      major_indexes.push_back((int)j);
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
    if(variants[i].gapped_freq >= 1-estimated_error && variants[i].cluster_assigned == -1){
      variants[i].consensus_numbers = all_genomes;
    }
  }

  //assign the number of the consensus genome
  for(uint32_t i=0; i < variants.size(); i++){
    for(uint32_t j=0; j < inverse_groups.size(); j++){
      //check to make sure you're lookin at a group that's part of the solution
      auto mit = std::find(solution.begin(), solution.end(), means[j]);
      if(mit == solution.end()) continue;

      //assign the point to all applicable groups
      auto it = std::find(inverse_groups[j].begin(), inverse_groups[j].end(), variants[i].cluster_assigned);
      if(it != inverse_groups[j].end()){
        variants[i].consensus_numbers.push_back(j);
      }
    }
  }
  amplicon_specific_cluster_assignment(variants, model);
  rewrite_position_masking(variants);
  std::vector<uint32_t> amplicons_to_mask;
  if(means.size() > 1){
    amplicons_to_mask = rewrite_amplicon_masking(variants, means);
  }
  modify_variant_masking(amplicons_to_mask, variants, means);
}
