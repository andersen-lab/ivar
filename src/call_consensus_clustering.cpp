#include "call_consensus_clustering.h"
#include "estimate_error.h"
#include "gmm.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>

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

bool account_for_clusters(std::vector<float> means, std::vector<std::vector<float>> results){
  bool keep = false;
  float error = 0.10;
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
    for(auto m : means){
      std::cerr << "remaining means " << m << std::endl;
    }
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

std::vector<std::vector<float>> find_solutions(std::vector<float> means){
  float error = 0.10;

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
    float sum = std::accumulate(results[i].begin(), results[i].end(), 0.0f);
    if(sum > 1-error && sum < 1+error){
      bool keep = true;
      if(keep){
        final_results.push_back(results[i]);
      }
    }
  }

  std::vector<float> useable_means;
  for(auto val : means){
    if(val != *min_iter && val != *max_iter){
      useable_means.push_back(val);
    }
  }
 
  //constrain that every solution must account for every cluster
  std::vector<std::vector<float>> final_final_results;
  
  for(uint32_t i=0; i < final_results.size(); i++){
    results.clear();
    current.clear();
    //our solution must contain more than one cluster (noise / universal min)
    if(final_results[i].size() < 2){
      continue;
    }
    //find combinations of the clusters
    find_combinations(final_results[i], 0, current, results);
    //account for points
    bool keep = account_for_clusters(useable_means, results);    
    if(keep){
      final_final_results.push_back(final_results[i]);
    }    
  }
    

  return(final_final_results);  
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

std::vector<std::vector<uint32_t>> find_combination_peaks(std::vector<float> solution, std::vector<float> means){ 
  std::vector<std::vector<uint32_t>> cluster_indexes(means.size());
  std::vector<float> current;
  std::vector<std::vector<float>> results;
  std::vector<float> totals;
  find_combinations(solution, 0, current, results);
  for(uint32_t i=0; i < results.size(); i++){
    float sum = std::accumulate(results[i].begin(), results[i].end(), 0.0f); 
    for(auto x : results[i]){
        std::cerr << x << " ";
    }
    std::cerr << "\n";
    std::cerr << "sum " << sum << std::endl;
    totals.push_back(sum);
  }
  //given a solution and the means, map each cluster to the cluster it contains
  for(uint32_t i=0; i < means.size(); i++){
    auto it = std::find(solution.begin(), solution.end(), means[i]);
    if(it != solution.end()){
        cluster_indexes[i].push_back(i);
    } else {
      float target = means[i];
      auto it = std::min_element(totals.begin(), totals.end(), [target](float a, float b) {return std::abs(a - target) < std::abs(b - target);});
      uint32_t index = std::distance(totals.begin(), it);
      for(auto x : results[index]){
        auto it2 = std::find(std::begin(means), std::end(means), x);
        uint32_t index2 = std::distance(std::begin(means), it2);
        cluster_indexes[i].push_back(index2);
      }
    }
  }
  /*for(uint32_t i=0; i < cluster_indexes.size(); i++){
    for(uint32_t j=0; j < cluster_indexes[i].size(); j++){
      std::cerr << cluster_indexes[i][j] << " ";
    }
    std::cerr << "\n";
  }*/
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
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, std::string variants_file){ 
  //output string
  float depth_cutoff = 10; 
  double error = 0.15; //acceptable error when forming solutions 
  
  float error_rate = cluster_error(variants_file);
  float freq_lower_bound = error_rate;
  float freq_upper_bound = 1 - error_rate;

  //read in the cluster values
  std::vector<float> means = parse_clustering_results(clustering_file);
  for(auto m : means){
    std::cerr << "consensus means " << m << std::endl;
  }
  
  std::vector<std::vector<float>> solutions = find_solutions(means);
  std::vector<float> solution;
  if(solutions.size() == 0){
    std::cerr << "no solution found" << std::endl;
    exit(1);
  } else if(solutions.size() > 1){
    std::cerr << "too many solutions" << std::endl;
    bool all_same = cluster_gravity_analysis(solutions);
    for(auto solution : solutions){
      for(auto s : solution){
        std::cerr << s << " ";
      }
      std::cerr << "\n";
    }
    if(all_same){
      solution = solutions[0];
    } else {
      exit(1);
    }
  } else{
    solution = solutions[0];
  }
  for(auto x : solution){
    std::cerr << x << std::endl;
  }
  std::vector<std::vector<uint32_t>> cluster_groups = find_combination_peaks(solution, means);
  std::vector<std::vector<uint32_t>> inverse_groups(means.size());
  for(uint32_t i=0; i < cluster_groups.size(); i++){
    for(uint32_t j=0; j < cluster_groups[i].size(); j++){
      //std::cerr << cluster_groups[i][j] << " i " << i << " j " << j << std::endl; 
      inverse_groups[cluster_groups[i][j]].push_back(i);
    }
  }
  /*for(auto x : inverse_groups){
    for(auto y : x){
        std::cerr << y << " ";
    }
    std::cerr << "\n";
  }
  exit(0);*/

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
  std::vector<int> minor_indexes;
  int max_index;
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
      if(tmp != largest){
        std::cerr << "adding to minor" << std::endl;
        max_index = (int)j;
        //minor_indexes.push_back((int)j);
      }        
      major_indexes.push_back((int)j);
    }
  }
  auto max_element = std::max_element(solution.begin(), solution.end());
  // Find the index of the largest element directly
  int index = std::distance(solution.begin(), max_element);
  float max_mean = solution[index];
  //if our largest cluster is under 0.50 we don't call
  if(max_mean < 0.50){
    std::cerr << "mean under 50" << std::endl;
    exit(1);
  }
  //find the largest position in the variants file
  uint32_t max_position = 0;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  } 
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
    if(variants[i].position == 1055){
      print = true;
      std::cerr << "\ntop freq " << variants[i].freq << " " << variants[i].nuc << " cluster " << variants[i].cluster_assigned << " " << variants[i].gapped_freq << std::endl;
      std::cerr << "vague assignment " << variants[i].vague_assignment << " del pos " << variants[i].pos_del_flag << " depth flag " << variants[i].depth_flag << std::endl;
      for(auto c : variants[i].probabilities){
        std::cerr << c << " ";
      }
      std::cerr << "\n";
    } else{
        print = false;
    }
    //if this position is experiencing fluctuation across amplicons, call ambiguity
    if(variants[i].amplicon_flux && variants[i].freq < freq_upper_bound){
      if(print){
        std::cerr << "position is experiencing fluctuation" << std::endl;
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
    std::cerr << i << " " << tmp[1054] << std::endl;
  }

  //write the consensus string to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);
  for(uint32_t i=0; i < all_sequences.size(); i++){
    float tmp_mean = means[i];
    auto it = std::find(solution.begin(), solution.end(), tmp_mean);
    if(it == solution.end()) continue;
    file << ">"+clustering_file+"_cluster_"+ std::to_string(means[i]) << "\n";
    file << all_sequences[i] << "\n";
  }
  file.close(); 
}
