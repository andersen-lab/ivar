#include "estimate_error.h"
#include "gmm.h"

std::vector<float> cluster_error(std::string filename){
  /*
    Here we use clustering to determine the value of the noise.
  */

  float lower_bound = 0.01;
  float upper_bound = 0.10;
  uint32_t depth_cutoff = 10;
  float quality_threshold = 20;
  uint32_t round_val = 4;
  
  std::vector<uint32_t> deletion_positions = find_deletion_positions(filename, depth_cutoff, lower_bound, upper_bound, round_val);
  std::vector<variant> base_variants;
  parse_internal_variants(filename, base_variants, depth_cutoff, lower_bound, upper_bound, deletion_positions, round_val);
  if(base_variants[0].version_1_var){
    calculate_reference_frequency(base_variants, filename, depth_cutoff, lower_bound, upper_bound, deletion_positions);
  }

  gaussian_mixture_model retrained_original;
  std::vector<variant> variants_original;
  uint32_t useful_count_original = 0;
  uint32_t max_pos = 0;
   for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].position > max_pos) max_pos = base_variants[i].position;
    //changed this for TEST
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].del_flag && !base_variants[i].amplicon_masked && !base_variants[i].primer_masked){
      
      useful_count_original++;
      variants_original.push_back(base_variants[i]);
    }
  }
  arma::mat data_original(1, useful_count_original, arma::fill::zeros);
  //std::cerr << useful_count_original << std::endl;
  uint32_t count_original=0;
  for(uint32_t i = 0; i < variants_original.size(); i++){
    //check if variant should be filtered for first pass model
    double tmp = static_cast<double>(variants_original[i].freq); //transform
    data_original.col(count_original) = tmp;
    count_original += 1;
  }
  uint32_t n = 3;
  retrained_original = retrain_model(n, data_original, variants_original, 2, 0.0001);
  assign_clusters(variants_original, retrained_original, 2);
  std::vector<std::vector<double>> clusters(n);
  for(auto var : variants_original){
    int cluster = var.cluster_assigned;
    clusters[cluster].push_back(var.freq);
    //std::cerr << var.freq << " " << cluster << std::endl;
  }
  uint32_t j = 0;
  uint32_t largest=0;
  std::vector<float> cluster_percents;
  for(uint32_t i=0; i < clusters.size(); i++){
    std::cerr << clusters[i].size() << std::endl;
    cluster_percents.push_back((float)clusters[i].size() / (float)variants_original.size());
    if(clusters[i].size() > largest){
      j = i;
      largest = clusters[i].size();
    }
  }
  for(auto me : retrained_original.means){
    std::cerr << me << std::endl;
  }  
  for(auto per : cluster_percents){
    std::cerr << per << std::endl;
  }
  //get the upper edge of the noise cluster
  auto max_it = std::max_element(clusters[j].begin(), clusters[j].end());
  std::cerr << *max_it << std::endl;
  std::vector<float> max_error_pos(max_pos+1);
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].freq <= *max_it){
      max_error_pos[base_variants[i].position] += base_variants[i].freq;
    }
  }

  float max_error = 0;  
  for(auto x : max_error_pos){
    if(x >= *max_it && x > max_error){
      max_error = x;
    }
  }
  //std::cerr << "max error " << max_error << std::endl;
  std::vector<float> tmp;
  //tmp.push_back(*max_it);
  //tmp.push_back(1-*max_it);

  tmp.push_back(max_error);
  tmp.push_back(1-max_error);
  //exit(0);
  return tmp;
}
