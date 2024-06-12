#include "call_consensus_clustering.h"
#include "gmm.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>

std::vector<double> pick_best_solution(std::vector<std::vector<double>> solutions){
  //which solution has the least error
  std::vector<double> errors;
  //which solution has the least number of populations
  std::vector<uint32_t> lengths;
  for(uint32_t i=0; i < solutions.size(); i++){
    double sum = std::accumulate(solutions[i].begin(), solutions[i].end(), 0.0);
    errors.push_back(sum);
    lengths.push_back(solutions[i].size());
  } 

  auto min_it = std::min_element(errors.begin(), errors.end());
  int min_index = std::distance(errors.begin(), min_it);
  return(solutions[min_index]);
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
  for(auto xx : means){
    std::cerr << xx << " ";
  }
  std::cerr << "\n";
 
  //TODO: Solve clustering
  std::vector<float> other_means;
  for(uint32_t m=0; m < means.size(); m++){
    std::cerr << "m " << means[m] << std::endl;
    if(means[m] < (float)0.97 && means[m] > (float)0.03){
      other_means.push_back(means[m]);
    }
  }
  std::vector<std::vector<double>> viable_solutions = solve_possible_solutions(other_means, error);
  std::vector<std::vector<double>> kept_solutions = deduplicate_solutions(viable_solutions);
  kept_solutions = deduplicate_solutions(viable_solutions);
  
  //we predicted a pure sample
  if(other_means.size() == 0){
    std::vector<double> tmp = {0.97, 0.03};
    kept_solutions.push_back(tmp);
  }
  if(kept_solutions.size() == 0){
    std::cerr << "no solutions found" << std::endl;
    exit(1);
  }
  for(auto xx : kept_solutions){
    std::cerr << "solution ";
    for(auto yy : xx){
      std::cerr << yy << " ";
    }
    std::cerr << "\n";
  }
  std::vector<double> solution;
  if(kept_solutions.size() > 1){
    //solution = pick_best_solution(kept_solutions);
    std::cerr << "too many solutions" << std::endl;
    exit(1);
  } else{ 
    solution = kept_solutions[0];
  }
  if(solution.size() > 3){
    exit(1);
  }
  std::vector<int> major_indexes;
  
  //index of the "100%" cluster
  for(uint32_t j=0; j < means.size(); j++){
    if(means[j] >= freq_upper_bound){
      major_indexes.push_back((int)j);
    }
  }
  auto max_it  = std::max_element(solution.begin(), solution.end());
  double largest_value = *max_it;
  auto it = std::find(means.begin(), means.end(), (float)largest_value);
  int index = static_cast<int>(std::distance(means.begin(), it));
  float max_mean = means[index];
  major_indexes.push_back(index);

  //if our largest cluster is under 0.50 we don't call
  if(max_mean < 0.50){
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
  //iterate all variants and determine
  for(uint32_t i = 0; i < variants.size(); i++){
   if(((float)variants[i].depth)*(1/variants[i].freq) < depth_cutoff){
     continue;
   } 
   if(variants[i].qual < 20 && variants[i].nuc != "-"){
     continue;
   }
   uint32_t position = variants[i].position;
   if(variants[i].low_prob_flag && means[variants[i].cluster_assigned] != (float)0.97 && variants[i].freq < max_mean){
      std::cerr << "low prob " << position << " " << variants[i].freq << " " << variants[i].nuc << " " << variants[i].probabilities[variants[i].cluster_assigned] << " " << variants[i].cluster_assigned << std::endl;
      std::cerr << means[variants[i].cluster_assigned] << std::endl;
      continue;
   }
   if(variants[i].vague_assignment && variants[i].freq < 0.97 && variants[i].freq < max_mean){      
     std::cerr << "vague assignment " << position << " " << variants[i].freq << " " << variants[i].nuc << std::endl;
    for(auto xx : variants[i].probabilities){
      std::cerr << xx << " ";
    }
    std::cerr << "\n";
    continue;
   }
   auto it = std::find(major_indexes.begin(), major_indexes.end(), variants[i].cluster_assigned);
   if(it != major_indexes.end()){
     consensus_vector[position-1] = variants[i].nuc;
    } else if(variants[i].freq > freq_upper_bound) {
      consensus_vector[position-1] = variants[i].nuc;
    }
  }
  //stitch to a consensus string 
  std::string consensus_sequence = std::accumulate(consensus_vector.begin(), consensus_vector.end(), std::string(""));
  //write the consensus string to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);
  file << ">"+clustering_file << "\n";
  file << consensus_sequence;
  file.close(); 
}
