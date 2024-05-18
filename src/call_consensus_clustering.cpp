#include "call_consensus_clustering.h"
#include "gmm.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>

int determine_assignment_status(variant var){
  if(!var.cluster_outlier && !var.low_prob_flag){
    return(var.cluster_assigned);
  }else{
    return(-1);
  }
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
  //techncially what we should do here is reload the model and reassign all variants including those that got excluded in the first pass
    
  //read in the cluster values
  std::vector<float> means = parse_clustering_results(clustering_file);
  /*for(auto xx : means){
    std::cerr << xx << " ";
  }
  std::cerr << "\n";*/
 
  //solve the clustering?
  //index of the "100%" cluster
  int universal_cluster = 0;
  for(uint32_t j=0; j < means.size(); j++){
    if(means[j] == freq_upper_bound){
      universal_cluster = (int)j;
    }
  }
  float largest_value = 0;
  int index_max_cluster = -2;
  //index the largest cluster
  for(uint32_t j=0; j < means.size(); j++) {
    if(means[j] > largest_value && means[j] != freq_upper_bound){
      largest_value = means[j];
      index_max_cluster = (int)j;
    }
  }
  if(means[index_max_cluster] < 0.50){
    index_max_cluster = -2;
  }

  //find the largest position in the variants file
  uint32_t max_position = 0;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  } 
  //std::cerr << index_max_cluster << " " << universal_cluster << std::endl;
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
   if(variants[i].vague_assignment){
     std::cerr << "vague assignment " << position << " " << variants[i].freq << " " << variants[i].nuc << " " << variants[i].amplicon_flux << " " << variants[i].primer_masked << std::endl;
    for(auto xx : variants[i].probabilities){
      std::cerr << xx << " ";
    }
    std::cerr << "\n";
   }
   if(variants[i].cluster_assigned == index_max_cluster){
      int assign = determine_assignment_status(variants[i]);
      if(assign == index_max_cluster){
        consensus_vector[position-1] = variants[i].nuc;
      }
    } else if(variants[i].cluster_assigned == universal_cluster || variants[i].freq > freq_upper_bound) {
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
