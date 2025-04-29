#include "estimate_error.h"
#include "call_consensus_clustering.h"
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

bool within_error_range(std::vector<float> values, float target, float error){
  //test if the sum of the vector equals the target value within some error
  float sum = std::accumulate(values.begin(), values.end(), 0.0f);
  if(sum < target+error && sum > target-error){
    return(true);  
  } else{
    return(false);
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
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means){ 
  std::cerr << "calling consensus" << std::endl;
  double max_mean=0;
  set_freq_range_flags(variants, 0, 1);
  double error_rate = cluster_error(variants, min_qual, min_depth);
  float freq_lower_bound = 1-error_rate+0.001;
  float freq_upper_bound = error_rate-0.001;
  set_freq_range_flags(variants, freq_lower_bound, freq_upper_bound);
  //find the largest position in the variants file
  uint32_t max_position = 0;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  }
  std::cerr << "lower " << freq_lower_bound << " upper " << freq_upper_bound << std::endl;
  bool print = false;
  uint32_t adjust_i = 0; //to help us track how many insertions we've seen
  uint32_t last_adjustment = 0;
  //initialize sequences for all possible populations
  std::vector<std::vector<std::string>> all_consensus_seqs;
  for(uint32_t i=0; i < means.size(); i++){
    std::vector<std::string> tmp(max_position, "N");
    all_consensus_seqs.push_back(tmp);
  }

  //order varaints by position
  std::sort(variants.begin(), variants.end(), [](const variant& a, const variant& b) {return a.position < b.position;}); 
  
  //iterate all variants and determine
  for(uint32_t i = 0; i < variants.size(); i++){
    //TESTLINES
    if(variants[i].position == 5){
      print = true;
      std::cerr << "\ntop freq " << variants[i].freq << " " << variants[i].nuc << " cluster " << variants[i].cluster_assigned << " " << variants[i].gapped_freq << std::endl;
      std::cerr << "vague assignment " << variants[i].vague_assignment << " depth flag " << variants[i].depth_flag << std::endl;
      std::cerr << "amplicon masked " << variants[i].amplicon_masked << " amp flux pos " << variants[i].amplicon_flux << std::endl;        
    }else{
      print = false;
    }
    float freq = variants[i].gapped_freq;
    float qual = variants[i].qual;
    uint32_t depth = variants[i].depth;    
    //depth, quality, and low frequency bypass
    if(freq < freq_lower_bound || qual < (float)min_qual || depth < min_depth){
      if(print) std::cerr << "min qual, freq, or depth issue " << qual << " " << freq << " " << depth << " flb " << freq_lower_bound << " mq " << (float)min_qual << " md " << min_depth << std::endl;
      continue;
    }

    //TODO handle the unresolved code portion

    //if this amplicon is experiencing fluctuation across amplicons, call ambiguity
    if(variants[i].amplicon_masked && variants[i].freq < freq_upper_bound){
      if(print){
        std::cerr << "amplicon is experiencing fluctuation" << std::endl;
      }
      continue;
    }
    //this variant position experiences fluctuation across amplicons
    if(variants[i].amplicon_flux){
      if(print){
        std::cerr << "amplicon in flux" << std::endl;
      }
      continue;
    } 
     uint32_t position = variants[i].position;
     if(variants[i].low_prob_flag && variants[i].freq < max_mean){
        if(print){
            std::cerr << "c" << std::endl;
            for(auto p : variants[i].probabilities){
                std::cerr << p << " ";
            }
            std::cerr << "\n";
        }
       continue;
     }
     if(variants[i].vague_assignment && variants[i].freq < freq_upper_bound && variants[i].freq < max_mean){      
       if(print){
          std::cerr << "d" << std::endl;
          for(auto a : variants[i].probabilities){
            std::cerr << a << " ";
          }
          std::cerr << "\n";
       }
       continue;
     }
     //handle all the cases where you never assigned anything, assign to all if it's over the upper bound
     if(variants[i].cluster_assigned == -1){
      if(variants[i].gapped_freq < freq_upper_bound) continue;
      if(print) std::cerr << "not assigned anything" << std::endl;
      for(uint32_t j=0; j < all_consensus_seqs.size(); j++){
        //TODO multiple insertions are not yet handled
        if(variants[i].nuc.find('+') != std::string::npos){
          std::string nuc = variants[i].nuc;
          nuc.erase(std::remove(nuc.begin(), nuc.end(), '+'), nuc.end());
          ++adjust_i;
          last_adjustment = variants[i].position;
          all_consensus_seqs[j].push_back("N"); //add an extra character
          all_consensus_seqs[j][position-1+adjust_i] = nuc;
        } else if (variants[i].position == last_adjustment){
          all_consensus_seqs[j][position-1+adjust_i-1] = variants[i].nuc;         
        } else {
          all_consensus_seqs[j][position-1+adjust_i] = variants[i].nuc;
        }
      }
      continue;
    }
    //TODO INVERSE CODE
    for(uint32_t j=0; j < variants[i].consensus_numbers.size(); j++){
        if(variants[i].nuc.find('+') != std::string::npos){
          std::string nuc = variants[i].nuc;
          ++adjust_i;
          last_adjustment = variants[i].position;
          all_consensus_seqs[variants[i].consensus_numbers[j]].push_back("N"); //add an extra character
          nuc.erase(std::remove(nuc.begin(), nuc.end(), '+'), nuc.end());
          all_consensus_seqs[variants[i].consensus_numbers[j]][position-1+adjust_i] = nuc;
          std::cerr << "position " << position-1+adjust_i << " " << adjust_i << std::endl;
        } else if (variants[i].position == last_adjustment){
          all_consensus_seqs[variants[i].consensus_numbers[j]][position-1+adjust_i-1] = variants[i].nuc;         
        } else {
          all_consensus_seqs[variants[i].consensus_numbers[j]][position-1+adjust_i] = variants[i].nuc;
        }
    }
  }
  std::vector<std::string> all_sequences;
  for(uint32_t i=0; i < all_consensus_seqs.size(); i++){
    std::string tmp = std::accumulate(all_consensus_seqs[i].begin(), all_consensus_seqs[i].end(), std::string(""));
    tmp.erase(std::remove(tmp.begin(), tmp.end(), '-'), tmp.end());
    all_sequences.push_back(tmp);
  }

  //write the consensus string to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);


  std::vector<uint32_t> indices(all_sequences.size());
  for (uint32_t i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }

  // Sort indices based on float values (descending)
  std::sort(indices.begin(), indices.end(), [&](uint32_t i, uint32_t j) {return means[i] > means[j];});

  // Apply sorted order
  std::vector<std::string> sorted_strings;
  std::vector<double> sorted_values;
  for (auto i : indices) {
    sorted_strings.push_back(all_sequences[i]);
    sorted_values.push_back(means[i]);
  }

  for(uint32_t i=0; i < sorted_strings.size(); i++){ 
    double tmp_mean = sorted_values[i];
    auto it = std::find(solution.begin(), solution.end(), tmp_mean);
    if(it == solution.end()){ 
      continue;
    }
    file << ">"+clustering_file+"_cluster_"+ std::to_string(tmp_mean) << "\n";
    file << sorted_strings[i] << "\n";
  }
  file.close(); 
}
