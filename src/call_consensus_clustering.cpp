#include "call_consensus_clustering.h"
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

  std::vector<std::vector<float>> final_results;
  //constrain that the solutions must add to 1
  for(uint32_t i=0; i < results.size(); i++){
    float sum = std::accumulate(results[i].begin(), results[i].end(), 0.0f);
    if(sum > 1-error && sum < 1+error){
      bool keep = true;
      //for(auto val : results[i]){
      //  if(val < 0.03 || val > 0.97) keep = false;
      //}
      if(keep){
        final_results.push_back(results[i]);
      }
    }
  }
  std::vector<float> useable_means;
  for(auto val : means){
    if(val > error || val > (1-error)){
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
void cluster_consensus(std::vector<variant> variants, std::string clustering_file){ 
  //output string
  float depth_cutoff = 10; 
  float freq_upper_bound = 0.97; 
  double error = 0.10; //acceptable error when forming solutions 
   
  //read in the cluster values
  std::vector<float> means = parse_clustering_results(clustering_file);
  for(auto m : means){
    std::cerr << m << std::endl;
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
  //populate a consensus vector with empty strings
  std::vector<std::string> consensus_vector(max_position, "N");  

  bool print = false;
  std::vector<uint32_t> skip_positions;
  //iterate all variants and determine
  for(uint32_t i = 0; i < variants.size(); i++){
    //auto it_skip = std::find(skip_positions.begin(), skip_positions.end(), variants[i].position);
    //if(it_skip != skip_positions.end()){
    //  continue;
    //}
    if(variants[i].position == 210){
      print = true;
      std::cerr << "\ntop freq " << variants[i].freq << " " << variants[i].nuc << " cluster " << variants[i].cluster_assigned << " " << variants[i].gapped_freq << std::endl;
      std::cerr << "vague assignment " << variants[i].vague_assignment << " del pos " << variants[i].pos_del_flag << std::endl;
      for(auto c : variants[i].probabilities){
        std::cerr << c << " ";
      }
      std::cerr << "\n";
    } else{
        print = false;
    }
    if(((float)variants[i].depth)*(1/variants[i].freq) < depth_cutoff){
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
     if(variants[i].low_prob_flag && means[variants[i].cluster_assigned] != (float)0.97 && variants[i].freq < max_mean){
        if(print){
            std::cerr << "c" << std::endl;
            for(auto p : variants[i].probabilities){
                std::cerr << p << " ";
            }
            std::cerr << "\n";
        }
       //std::cerr << "low prob " << position << " " << variants[i].freq << " " << variants[i].nuc << " " << variants[i].probabilities[variants[i].cluster_assigned] << " " << variants[i].cluster_assigned << std::endl;
       //std::cerr << means[variants[i].cluster_assigned] << std::endl;
       continue;
     }
     if(variants[i].vague_assignment && variants[i].freq < 0.97 && variants[i].freq < max_mean){      
       if(print){
          std::cerr << "d" << std::endl;
          //std::cerr << "vague assignment " << position << " " << variants[i].freq << " " << variants[i].nuc << std::endl;
          for(auto a : variants[i].probabilities){
            std::cerr << a << " ";
          }
          std::cerr << "\n";
       }
       continue;
     }
     auto it = std::find(major_indexes.begin(), major_indexes.end(), variants[i].cluster_assigned);
     if(it != major_indexes.end()){
      if(print){
         std::cerr << "in major cluster" << std::endl;
       }
       consensus_vector[position-1] = variants[i].nuc;
       //sometimes we assign a nuc and a deletion to consensus, due to the effect of ungapped depth on frequency, however if we assign a deletion we should not try to assign anything else
       /*if(variants[i].nuc == "-"){
         skip_positions.push_back(position);
       }*/
     } else if(variants[i].pos_del_flag && variants[i].gapped_freq > freq_upper_bound) {
       if(print){
         std::cerr << "del flag greater than upper bound " << position << " " << variants[i].nuc << std::endl;
       }
       consensus_vector[position-1] = variants[i].nuc;
     } else if(!variants[i].pos_del_flag && variants[i].freq > freq_upper_bound){
       if(print){
         std::cerr << "greated than upper bound " << position << " " << variants[i].nuc << std::endl;
       }
       consensus_vector[position-1] = variants[i].nuc;
     }
  }
  //stitch to a consensus string 
  std::string consensus_sequence = std::accumulate(consensus_vector.begin(), consensus_vector.end(), std::string(""));

  //write the consensus string to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);
  file << ">"+clustering_file << "\n";
  file << consensus_sequence << "\n";
  file.close(); 
}
