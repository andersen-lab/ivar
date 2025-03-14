#include "call_consensus_clustering.h"
#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>

bool test_cluster_deviation(float nearest_cluster, float variant_cluster, float std_dev){
  bool fluctuation = false;
  //CLEANUP THIS CAN BE CALCULATED ONCE PER ALL CLUSTERS
  //determine if the assigned and nearest cluster can be resolved based on variant fluctuation
  std::vector<double> tmp = {(double) nearest_cluster, (double) variant_cluster};
  double cluster_dev = calculate_standard_deviation(tmp);
  if((double)std_dev > cluster_dev){
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

void call_majority_consensus(std::vector<variant> variants, uint32_t max_position, std::string clustering_file, double default_threshold){
  //if we can't find a solution simply take the majority variant per position
  std::vector<std::string> nucs;
  std::vector<float> freqs;
  std::vector<std::string> tmp(max_position, "N");
  for(uint32_t i=1; i <= max_position; i++){
    freqs.clear();
    nucs.clear();
    for(uint32_t j=0; j < variants.size(); j++){
      if(variants[j].position == i){
        nucs.push_back(variants[j].nuc);
        freqs.push_back(variants[j].freq);
      }
    }
    if(freqs.size() == 0) continue;
    uint32_t index = std::distance(freqs.begin(), std::max_element(freqs.begin(), freqs.end()));
    if(freqs[index] >= (float)default_threshold){
      tmp[i-1] = nucs[index];
    }
  }
  std::string consensus_string = std::accumulate(tmp.begin(), tmp.end(), std::string(""));
  //write the consensus to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);
  std::string name = ">"+clustering_file+"_"+std::to_string(default_threshold)+"_threshold";
  file << name << "\n";
  file << consensus_string << "\n";
  file.close(); 
}

float find_nearest_distance(const std::vector<float> all_sums, float value) {
    float min_distance = std::numeric_limits<float>::max();
    for (auto num : all_sums) {
        float distance = std::abs(num - value);
        if (distance < min_distance) {
            min_distance = distance;
        }
    }
    return min_distance;
}

bool account_peaks(std::vector<float> possible_solution, std::vector<float> means, float total, float error){
  bool valid = true;
  std::vector<float> current;
  std::vector<std::vector<float>> results;
  find_combinations(possible_solution, 0, current, results);
 
  std::vector<float> all_sums; 
  for(auto result : results){
    float sum = std::accumulate(result.begin(), result.end(), 0.0f);
    all_sums.push_back(sum);
  }

  //check if all means can be accounted for
  for(auto mean : means){
    float dist = find_nearest_distance(all_sums, mean);
    if(dist > error){
      valid = false;
      break;
    }
  }
  return(valid);
}

bool within_error_range(std::vector<float> values, float target, float error){
  //test if the sum of the vector equals the target value within some error
  float sum = std::accumulate(values.begin(), values.end(), 0.0f);
  if(sum < target+error && sum > target-error){
    return(true);  
  } else{
    return(false);
  }
}

std::vector<std::vector<float>> find_subsets_with_error(std::vector<float> means, float target, float error){
  //first we find all the possible combinations
  std::vector<float> current;
  std::vector<std::vector<float>> results;
  find_combinations(means, 0, current, results);

  std::vector<std::vector<float>> valid_combinations;  
  for(uint32_t i=0; i < results.size(); i++){
    bool in_range = within_error_range(results[i], target, error);
    if(in_range){
      valid_combinations.push_back(results[i]);
    }
  }
  return(valid_combinations);
}

std::vector<std::vector<float>> frequency_pair_finder(std::vector<variant> variants, float lower_bound, float upper_bound, std::vector<float> means){ 
  std::vector<std::vector<float>> pairs;
  std::vector<uint32_t> track_positions;

  for(uint32_t i=0; i < variants.size(); i++){
    if(!variants[i].depth_flag && !variants[i].qual_flag && !variants[i].outside_freq_range && variants[i].cluster_assigned != -1){
      auto it = std::find(track_positions.begin(), track_positions.end(), variants[i].position);
      //found
      if(it != track_positions.end()){
        size_t index = std::distance(track_positions.begin(), it);
        pairs[index].push_back(means[variants[i].cluster_assigned]);
      } else{
        std::vector<float> tmp = {means[variants[i].cluster_assigned]};
        pairs.push_back(tmp);
        track_positions.push_back(variants[i].position);
      }
    }    
  } 

  return(pairs);
}

bool cluster_gravity_analysis(std::vector<std::vector<float>> solutions){
  //in the event of multiple solutions, check that the largest cluster is the same
  std::vector<float> max_values;
  for(auto solution : solutions){
    float max = *std::max_element(solution.begin(), solution.end());
    max_values.push_back(max);
  }
  bool all_same = std::all_of(max_values.begin() + 1, max_values.end(), [&](float x) { return x == max_values[0]; });  
  return(all_same);
}

bool account_for_clusters(std::vector<float> means, std::vector<std::vector<float>> results, float error){
  bool keep = false;
  std::vector<float> accounted_means;

  for(uint32_t i=0; i < results.size(); i++){
    float total = std::accumulate(results[i].begin(), results[i].end(), 0.0f);
    //determine if this is close to a cluster
    for(uint32_t j=0; j < means.size(); j++){
      float diff = std::abs(total-means[j]);
      if(diff < error){
        accounted_means.push_back(means[j]);
      }
    }
  }
  
  for(auto val : accounted_means){
    auto it = std::find(means.begin(), means.end(), val);
    if (it != means.end()){
      uint32_t index = std::distance(means.begin(), it);
      means.erase(means.begin() + index);
    }
  }  
  if(means.size() == 0){
    keep = true;
  } else{
    keep = false;
  }
  return(keep);
}

void find_combinations(std::vector<float> means, uint32_t index, std::vector<float> &current, std::vector<std::vector<float>> &results){
  if (!current.empty()){
    results.push_back(current);
  }
  for (uint32_t i = index; i < means.size(); ++i) {
    current.push_back(means[i]);
    find_combinations(means, i+1, current, results);
    current.pop_back();
  }
}

std::vector<std::vector<float>> find_solutions(std::vector<float> means, float error){
  std::vector<float> current;
  std::vector<std::vector<float>> results;
  find_combinations(means, 0, current, results);
  
  std::sort(results.begin(), results.end());
  results.erase(std::unique(results.begin(), results.end()), results.end());

  auto max_iter = std::max_element(means.begin(), means.end());
  auto min_iter = std::min_element(means.begin(), means.end());

  std::vector<std::vector<float>> final_results;
  //constrain that the solutions must add to 1
  for(uint32_t i=0; i < results.size(); i++){
    bool keep = within_error_range(results[i], 1, error);
    if(keep){
      final_results.push_back(results[i]);
    }
  }
  return(final_results);  
}

std::vector<float> parse_string_to_vector(const std::string& str) {
    std::vector<float> result;
    std::stringstream ss(str);
    char ch; // Used to read and discard non-numeric characters, including the decimal point
    float num;

    // Read characters one by one
    while (ss >> ch) {
        // Check if the character is a digit, a minus sign, or a decimal point
        if ((ch >= '0' && ch <= '9') || ch == '-' || ch == '.') {
            // Put back the character into the stream to correctly read the number
            ss.putback(ch);
            ss >> num; // Read the number as float
            result.push_back(num); // Add the number to the vector
        }
    }

    return result;
}

std::vector<std::vector<uint32_t>> find_combination_peaks(std::vector<float> solution, std::vector<float> means, std::vector<float> &unresolved){ 
  std::vector<std::vector<uint32_t>> cluster_indexes(means.size());
  std::vector<float> current;
  std::vector<std::vector<float>> results;
  std::vector<float> totals;

  find_combinations(solution, 0, current, results);
  for(uint32_t i=0; i < results.size(); i++){
    float sum = std::accumulate(results[i].begin(), results[i].end(), 0.0f); 
    totals.push_back(sum);
  }
  //given a solution and the means, map each cluster to the cluster it contains
  for(uint32_t i=0; i < means.size(); i++){
    auto it = std::find(solution.begin(), solution.end(), means[i]);
    //th mean is part of the solution
    if(it != solution.end()){
        cluster_indexes[i].push_back(i);
        float target = means[i];
        std::vector<float> distances(totals.size());
        std::transform(totals.begin(), totals.end(), distances.begin(), [target](float num) { return std::abs(target - num); }); 
        uint32_t count = 0;
        
        //this checks the distances from the mean to all other possible peaks
        for(uint32_t d=0; d < distances.size(); d++){
          if(distances[d] < 0.03) count += 1;
        }
        if(count > 1) unresolved.push_back(target);
       
    } else {
      float target = means[i];
      //the problem with this is that it looks at the min but not if two overlapping peaks occur
      auto it = std::min_element(totals.begin(), totals.end(), [target](float a, float b) {return std::abs(a - target) < std::abs(b - target);});
      
      std::vector<float> distances(totals.size());
      std::transform(totals.begin(), totals.end(), distances.begin(), [target](float num) { return std::abs(target - num); }); 
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
  }
  for(auto u : unresolved) std::cerr << u << std::endl;*/
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

std::vector<float> parse_clustering_results(std::string clustering_file){
  std::ifstream infile(clustering_file + ".txt");
  std::string line;
  uint32_t count = 0;
  std::vector<float> numbers;
  while (std::getline(infile, line)) {
    if(count == 0) {
      count += 1;
      continue;
    }
    std::vector<std::string> row_values;
    split(line, '\t', row_values);
    std::string means = row_values[0];
    numbers = parse_string_to_vector(means);
    count += 1;
  }  
  return(numbers);
}
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, std::string variants_file, double default_threshold){ 
  float depth_cutoff = 10; 
  double error = 0.10; 
  float solution_error = 0.05;
  double quality_threshold = 20; 

  std::vector<float> error_rate = cluster_error(variants_file, quality_threshold);
  float freq_lower_bound = error_rate[0];
  float freq_upper_bound = error_rate[1];

  //read in the cluster values
  std::vector<float> means = parse_clustering_results(clustering_file);
  for(auto m : means){
    std::cerr << "consensus means " << m << std::endl;
  }
  std::vector<std::vector<float>> clusters(means.size());
  for(auto var : variants){
    if(var.cluster_assigned != -1){
      clusters[var.cluster_assigned].push_back(var.freq);
    }
  }
  //find the largest position in the variants file
  uint32_t max_position = 0;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  }
  //find position wise frequency pairs
  std::vector<std::vector<float>> pairs = frequency_pair_finder(variants, freq_lower_bound, freq_upper_bound, means); 
  std::vector<std::vector<float>> solutions = find_solutions(means, error);  
  
  //find peaks that can't be a subset of other peaks
  std::vector<float> non_subset_means;
  for(uint32_t i=0; i < means.size(); i++){
    std::vector<std::vector<float>> tmp = find_subsets_with_error(means, means[i], solution_error);    
    if(tmp.size() <= 1){
      non_subset_means.push_back(means[i]);
    }
  }
  //reduce solution space to things that contain the non subset peaks
  std::vector<std::vector<float>> realistic_solutions;
  for(uint32_t i=0; i < solutions.size(); i++){  
      std::vector<float> tmp = solutions[i];
      bool found = std::all_of(non_subset_means.begin(), non_subset_means.end(), [&tmp](float value) {return std::find(tmp.begin(), tmp.end(), value) != tmp.end();});
     if(found){
        realistic_solutions.push_back(solutions[i]);
     }
  }
  //check each solution that every possible peak is accounted for
  std::vector<std::vector<float>> solution_sets;
  for(uint32_t i=0; i < realistic_solutions.size(); i++){
    bool keep = account_peaks(realistic_solutions[i], means, 1, solution_error);
    if(keep){
      solution_sets.push_back(realistic_solutions[i]);
    }
  }

  for(auto sol : solution_sets){
    std::cerr << "\nsolution" << std::endl;
    for(auto s : sol){
      std::cerr << s << " ";
    }
  }
  std::cerr << "\n" << std::endl;

  std::vector<float> solution;
  bool traditional_majority= false; //if we can't find a solution call a traditional majority consensus
  if(solution_sets.size() == 0){
    std::cerr << clustering_file << " no solution found" << std::endl;
    traditional_majority = true;
  } else if(solution_sets.size() > 1){
    std::cerr << "too many solutions" << std::endl;
    traditional_majority = true;
  } else{
    solution = solution_sets[0];
  }

  if(traditional_majority){
    call_majority_consensus(variants, max_position, clustering_file, default_threshold);
    exit(0);
  }

  for(auto x : solution){
    std::cerr << x << std::endl;
  }
  std::vector<float> unresolved;
  std::vector<std::vector<uint32_t>> cluster_groups = find_combination_peaks(solution, means, unresolved);
  
  std::vector<std::vector<uint32_t>> inverse_groups(means.size());
  for(uint32_t i=0; i < cluster_groups.size(); i++){
    for(uint32_t j=0; j < cluster_groups[i].size(); j++){
      inverse_groups[cluster_groups[i][j]].push_back(i);
    }
  }
  //TESTLINES MEAN CODE
  std::string solution_string = "[";
  for(uint32_t i=0; i < solution.size(); i++){
    if(i != 0){
      solution_string += ",";
    }
    std::string tmp = std::to_string(solution[i]);
    solution_string += tmp;
  }

  solution_string += "]";
  std::string solution_filename = clustering_file + "_solution.txt";
  std::ofstream file_sol(solution_filename);
  file_sol << "means\n";
  file_sol << solution_string << "\n";
  file_sol.close();

  float largest = *std::max_element(solution.begin(), solution.end());
  //define the clusters which contain the majority population
  std::vector<std::vector<float>> possible_clusters;
  std::vector<float> current;
  find_combinations(solution, 0, current, possible_clusters); 
  std::vector<float> expected_clusters; 
  for(uint32_t i=0; i < possible_clusters.size(); i++){
    bool keep = false;
    for(uint32_t j=0; j < possible_clusters[i].size(); j++){
      if(possible_clusters[i][j] == largest){
        keep = true;
        break;
      }
    }
    if(keep){
      float sum = std::accumulate(possible_clusters[i].begin(), possible_clusters[i].end(), 0.0f);
      expected_clusters.push_back(sum);
    }
  } 
  //a list of cluster assignments that we assign to consensus
  std::vector<int> major_indexes;
  //index of the "100%" cluster
  for(uint32_t j=0; j < means.size(); j++){
    float tmp = means[j];
    auto closest = *std::min_element(expected_clusters.begin(), expected_clusters.end(), [tmp](float a, float b) {
      return std::abs(a - tmp) < std::abs(b - tmp);
    });
    float diff = std::abs(closest - tmp); 
    auto it = std::find(solution.begin(), solution.end(), tmp);
    if((diff < error && it == solution.end()) || tmp == largest){
      std::cerr << "major index " << means[j] << " " << j << std::endl;
      major_indexes.push_back((int)j);
    }
  }
  auto max_element = std::max_element(solution.begin(), solution.end());
  // Find the index of the largest element directly
  int index = std::distance(solution.begin(), max_element);
  float max_mean = solution[index];

  bool print = false;
  std::vector<std::vector<std::string>> all_consensus_seqs;
  for(uint32_t i=0; i < means.size(); i++){
    std::vector<std::string> tmp(max_position, "N");
    all_consensus_seqs.push_back(tmp);
  }
   
  //iterate all variants and determine
  for(uint32_t i = 0; i < variants.size(); i++){
    //TESTLINES
    if(variants[i].nuc.find('+') != std::string::npos) continue;
    //TESTLINES
    if(variants[i].position == 13572){
      print = true;
      std::cerr << "\ntop freq " << variants[i].freq << " " << variants[i].nuc << " cluster " << variants[i].cluster_assigned << " " << variants[i].gapped_freq << std::endl;
      std::cerr << "vague assignment " << variants[i].vague_assignment << " del pos " << variants[i].pos_del_flag << " depth flag " << variants[i].depth_flag << std::endl;
      std::cerr << variants[i].amplicon_masked << std::endl;
    } else{
        print = false;
    }
    //if the mean for this cluster is unresolved we skip it
    auto it = std::find(unresolved.begin(), unresolved.end(), means[variants[i].cluster_assigned]);      
    if(it != unresolved.end()){ 
      if(print){
        std::cerr << "unresolved " << means[variants[i].cluster_assigned] << std::endl;
      }
      continue;
    }
    //if this position is experiencing fluctuation across amplicons, call ambiguity
    if(variants[i].amplicon_flux && variants[i].freq < freq_upper_bound && variants[i].freq > freq_lower_bound){

      //find all clusters not part of the same assignment
      std::vector<double> other_population_clusters;
      for(uint32_t j=0; j < inverse_groups.size(); j++){
        //check to make sure you're lookin at a group that's part of the solution
        auto mit = std::find(solution.begin(), solution.end(), means[j]);
        if(mit == solution.end()) continue;
        auto it = std::find(inverse_groups[j].begin(), inverse_groups[j].end(), variants[i].cluster_assigned);      
        //assigned cluster is not apart of the population
        if(it == inverse_groups[j].end())
        for(auto ig : inverse_groups[j]){
          //CLEAN UP this will push redundant things back
          other_population_clusters.push_back(means[ig]);
        }
      }

      //find the second closest cluster index
      double closest_mean = find_neighboring_cluster((double)variants[i].gapped_freq, variants[i].cluster_assigned, other_population_clusters);
      //check if the cluster is within the standard dev of the variant
      bool fluctuating = test_cluster_deviation(closest_mean, means[variants[i].cluster_assigned], variants[i].std_dev);       
      if(print){
        std::cerr << "fluctuating " << fluctuating << std::endl;
      }
      if(fluctuating){
        continue;
      }
    }

    //if this amplicon is experiencing fluctuation across amplicons, call ambiguity
    if(variants[i].amplicon_masked && variants[i].freq < freq_upper_bound){
      if(print){
        std::cerr << "amplicon is experiencing fluctuation" << std::endl;
      }
      continue;
    }

    if(variants[i].depth_flag){
      if(print){
        std::cerr << "a " << variants[i].depth_flag << std::endl;
      }
      continue;
     } 
     if(variants[i].qual < 20 && variants[i].nuc != "-"){
        if(print){
          std::cerr << "b" << std::endl;
       }
       continue;
     }
     uint32_t position = variants[i].position;
     if(variants[i].low_prob_flag && means[variants[i].cluster_assigned] != freq_upper_bound && variants[i].freq < max_mean){
        if(print){
            std::cerr << "c" << std::endl;
            for(auto p : variants[i].probabilities){
                std::cerr << p << " ";
            }
            std::cerr << "\n";
        }
       continue;
     }
     if(variants[i].vague_assignment && variants[i].freq < freq_upper_bound && variants[i].freq < max_mean && std::abs(variants[i].freq - means[variants[i].cluster_assigned]) > 0.10){      
       if(print){
          std::cerr << "d" << std::endl;
          for(auto a : variants[i].probabilities){
            std::cerr << a << " ";
          }
          std::cerr << "\n";
       }
       continue;
     }
     if(variants[i].freq <= freq_lower_bound) continue;
     //handle all the cases where you never assigned anything
    if(variants[i].cluster_assigned == -1){
      if(variants[i].pos_del_flag && variants[i].gapped_freq < freq_upper_bound) continue;
      if(!variants[i].pos_del_flag && variants[i].freq < freq_upper_bound) continue;
      //ADD IN ADDING TO ALL CLUSTERS
      for(uint32_t j=0; j < all_consensus_seqs.size(); j++){
        all_consensus_seqs[j][position-1] = variants[i].nuc;
      }
      continue;
    }

    for(uint32_t j=0; j < inverse_groups.size(); j++){
      //check to make sure you're lookin at a group that's part of the solution
      auto mit = std::find(solution.begin(), solution.end(), means[j]);
      if(mit == solution.end()) continue;
      //assign the point to all applicable groups
      auto it = std::find(inverse_groups[j].begin(), inverse_groups[j].end(), variants[i].cluster_assigned);      
      if(it != inverse_groups[j].end()){ 
        all_consensus_seqs[j][position-1] = variants[i].nuc;
      }
    }
  }
  std::vector<std::string> all_sequences;
  for(uint32_t i=0; i < all_consensus_seqs.size(); i++){
    std::string tmp = std::accumulate(all_consensus_seqs[i].begin(), all_consensus_seqs[i].end(), std::string(""));
    all_sequences.push_back(tmp);
  }

  //write the consensus string to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);
  for(uint32_t i=0; i < all_sequences.size(); i++){
    float tmp_mean = means[i];
    auto it = std::find(solution.begin(), solution.end(), tmp_mean);
    if(it == solution.end()) continue;
    file << ">"+clustering_file+"_cluster_"+ std::to_string(means[i]) << "\n";
    std::cerr << all_sequences[i][13571] << std::endl;
    file << all_sequences[i] << "\n";
  }
  file.close(); 
}
