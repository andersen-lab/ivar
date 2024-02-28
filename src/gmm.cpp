#include "./include/armadillo"
#include "gmm.h"
#include "call_consensus_clustering.h"
#include "population_estimate.h"
#include <fstream>
#include <cmath>
#include <algorithm>

bool is_substring(const std::string& main_string, const std::string& sub_string) {
    // Find returns std::string::npos if the substring is not found
    return main_string.find(sub_string) != std::string::npos;
}

double calculate_median(const std::vector<double>& data) {
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());

    size_t size = sorted_data.size();

    if (size % 2 == 0) {
        // If the size is even, return the average of the middle two elements
        return (sorted_data[size / 2 - 1] + sorted_data[size / 2]) / 2.0;
    } else {
        // If the size is odd, return the middle element
        return sorted_data[size / 2];
    }
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

void generate_combinations(std::vector<double> &input, std::vector<double>& current_combination, uint32_t start_index, uint32_t length, std::vector<std::vector<double>> &collect_combos) {
  if (length == 0) {
    collect_combos.push_back(current_combination);
    return;
  }

  for (uint32_t i = start_index; i < input.size(); i++) {
    current_combination.push_back(input[i]);
    generate_combinations(input, current_combination, i + 1, length - 1, collect_combos);
    current_combination.pop_back();
  }
}


std::vector<std::vector<double>> remove_unexplainable_solutions(std::vector<std::vector<double>> solutions, std::vector<double> means){
  double error = 0.10;
  std::vector<std::vector<double>> kept_solutions;
  for(auto vec : solutions){
    std::vector<double> possible_peaks;
    for(auto x : vec){
      //std::cerr << x << " ";
      possible_peaks.push_back(x);
    }
    //std::cerr << std::endl;
    std::vector<std::vector<double>> all_combinations; 
    for (uint32_t i= 1; i <= vec.size(); i++) {
      std::vector<double> current_combination;
      generate_combinations(vec, current_combination, 0, i, all_combinations);  
    }
    for(auto vec_two : all_combinations){
      double sum = std::accumulate(vec_two.begin(), vec_two.end(), 0.0); 
      possible_peaks.push_back(sum);
    }
    
    //find each mean
    bool keep_solution = true;
    for(auto mean : means){
      bool useful = false;
      if(mean < error){
        continue;
      }
      for(auto theoretical_peak : possible_peaks){
        //std::cerr << "theoretical peak " << theoretical_peak << " mean " << mean << std::endl;
        if(std::abs(theoretical_peak - mean) < error){
          //std::cerr << "this " << theoretical_peak << " " << mean << std::endl; 
          useful = true;
          break;
        }  
      }
      if(!useful){
        //std::cerr << mean << std::endl;
        keep_solution = false;
        break;          
      }
    }
    if(keep_solution){
      kept_solutions.push_back(vec);
    }
  }
  return(kept_solutions);
}

std::vector<std::vector<double>> solve_possible_solutions(std::vector<double> tmp_means){
  uint32_t n = 6;
  double error = 0.10;
  std::vector<double> means;
  for(auto &x : tmp_means){
    for(uint32_t i=1; i < n; i++){
      if(x * i < 1+error){
        means.push_back(x);
      }
    }
  }
  /*for(auto x : means){
    std::cerr << x << std::endl;
  }
  exit(1);*/
  std::vector<std::vector<double>> all_combinations;
  for (uint32_t i= 1; i < means.size(); i++) {
    std::vector<double> current_combination;
    generate_combinations(means, current_combination, 0, i, all_combinations);
    
  }
  std::vector<std::vector<double>> viable_solutions;
  for(auto vec : all_combinations){
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0); 
    /*for(auto x : vec){
      std::cerr << x << " ";
    }
    std::cerr << "sum " << sum << std::endl;*/
    if(sum < 1+error && sum > 1-error && vec.size() > 1) {
      //std::cerr << "saved" << std::endl;
      viable_solutions.push_back(vec);
    }
  }
  //check that these solutions explain all peaks
  std::vector<std::vector<double>> kept_solutions = remove_unexplainable_solutions(viable_solutions, means);  
  return(kept_solutions);
}

std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>>& input_vector) {
  std::vector<std::vector<double>> transposed_vector;
  // Check if the input vector is not empty
  if (!input_vector.empty() && !input_vector[0].empty()) {
    size_t rows = input_vector.size();
    size_t cols = input_vector[0].size();
    // Resize the transposed vector
    transposed_vector.resize(cols, std::vector<double>(rows));

    // Transpose the matrix
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        transposed_vector[j][i] = input_vector[i][j];
      }
    }
  }
  return transposed_vector;
}

double calculate_cluster_bounds(std::vector<variant> variants, uint32_t n){
  /*
   * Find the smallest point in the largest cluster and determine the threshold to be just below.
   */
  std::vector<std::vector<double>> clusters;
  for(uint32_t i=0; i < n; i++){
    std::vector<double> tmp;
    clusters.push_back(tmp);
  }
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].cluster_assigned > -1){
      clusters[variants[i].cluster_assigned].push_back(variants[i].freq);
    }
  }
  uint32_t max_idx = 0;
  double max_mean;
  for(uint32_t i=0; i < clusters.size(); i++){
    double sum = std::accumulate(clusters[i].begin(), clusters[i].end(), 0.0);
    double mean = sum / clusters[i].size();
    if(mean > max_mean){
      max_idx = i;
      max_mean = mean;
    }
  }
  double min_freq=1.0;
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].cluster_assigned == (int)max_idx && variants[i].freq < min_freq){
      if(!variants[i].cluster_outlier){
        min_freq = variants[i].freq;
      }
    }
  }
  return(min_freq);  
}

uint32_t count_useful_variants(std::vector<variant> variants){
  uint32_t count = 0;
  //determine the number of variants useful for modeling
  for(uint32_t i=0; i< variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
      count += 1;
    }
  }
  return(count);
}

void perm_generator(int n, int k, std::vector<std::vector<uint32_t>> &possible_permutations){
    std::vector<uint32_t> d(n);
    std::iota(d.begin(),d.end(),0);
    do {
        std::vector<uint32_t> tmp;
        for (int i = 0; i < k; i++){
          tmp.push_back(d[i]);
        }
        possible_permutations.push_back(tmp);
        std::reverse(d.begin()+k,d.end());
    } while(std::next_permutation(d.begin(),d.end()));
}

std::vector<uint32_t> compare_cluster_assignment(std::vector<std::vector<double>> prob_matrix, std::vector<uint32_t> assigned){
  double threshold = 3;
  std::vector<uint32_t> flagged_idx;
  for(uint32_t i=0; i < prob_matrix.size(); i++){
    double assigned_prob = prob_matrix[i][assigned[i]];
    std::vector<double> tmp = prob_matrix[i];
    std::sort(tmp.begin(), tmp.end(), std::greater<double>());
    for(uint32_t j=0; j < tmp.size(); j++){
      if(tmp[j] >= assigned_prob) continue;
        if(exp(tmp[j]) * threshold > exp(assigned_prob)){
          //std::cerr << exp(tmp[j]) << " " << exp(assigned_prob) << std::endl;
          flagged_idx.push_back(i);
        }
        break;
    }
  }
  return(flagged_idx);
}

std::vector<uint32_t> calculate_joint_probabilities(std::vector<std::vector<double>> prob_matrix, std::vector<std::vector<uint32_t>> permutations){
  /*
   * Calcualte the best assignments maximizing the probability. Prob matrix is (n_clusters, n_variants)
   */
  std::vector<double> scores; //score for every permuation of assignment
  //on the permutation level                            
  for(uint32_t i=0; i < permutations.size(); i++){   
   double score = 0; 
    for(uint32_t j = 0; j < permutations[i].size(); j++){
      //num of variants in position must match permutation
      if(permutations[i].size() != prob_matrix.size()){
        score = -1000000;
        continue;
      }
      score += prob_matrix[j][permutations[i][j]];
    }
    scores.push_back(score);
  }
  uint32_t max_idx = std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));
  return(permutations[max_idx]);
}

void assign_variants_simple(std::vector<variant> &variants, std::vector<std::vector<double>> prob_matrix){
  uint32_t n = prob_matrix.size();

  //find the unique positions
  std::vector<uint32_t> unique_pos;
  for(uint32_t i = 0; i < variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
      if (std::find(unique_pos.begin(), unique_pos.end(), variants[i].position) == unique_pos.end()) {
        unique_pos.push_back(variants[i].position);
      }
    }
  }
  //determine all possible permutations of assignments, disregard sum condition of E(u) ~= 1
  std::vector<uint32_t> idx_combinations(n);
  std::iota(idx_combinations.begin(), idx_combinations.end(), 0);
  std::vector<std::vector<uint32_t>> possible_permutations;

  perm_generator(n, 1, possible_permutations);
  perm_generator(n, 2, possible_permutations);
  perm_generator(n, 3, possible_permutations);
  perm_generator(n, 4, possible_permutations);
   
  //now we loop every unique position and assign the max prob combo of variants
  for(uint32_t i=0; i < unique_pos.size(); i++){
    std::vector<uint32_t> pos_idxs;
    std::vector<std::vector<double>> tmp_prob;
    uint32_t j = 0;
    //all locations in the prob matrix for this position
    for(uint32_t k = 0;  k < variants.size(); k++){
      if(variants[k].amplicon_flux || variants[k].depth_flag || variants[k].outside_freq_range || variants[k].qual_flag || variants[k].del_flag) continue;
      
      if(variants[k].position  == unique_pos[i]){
        pos_idxs.push_back(j);
        std::vector<double> tmp;
        for(uint32_t l=0; l < n; l++){
          tmp.push_back(prob_matrix[l][j]);
        }
        tmp_prob.push_back(tmp);
      }
      j++;
    }
    //assign variants based on most probable position-wise
    
    std::vector<uint32_t> assigned = calculate_joint_probabilities(tmp_prob, possible_permutations);
    //make sure the assignment is concrete
    std::vector<uint32_t> assignment_flagged = compare_cluster_assignment(tmp_prob, assigned);
    for(uint32_t j=0; j < pos_idxs.size(); j++){
      std::vector<uint32_t>::iterator tmp = std::find(assignment_flagged.begin(), assignment_flagged.end(), j);
      uint32_t k = 0;
      for(uint32_t z =0; z < variants.size(); z++){
        if(variants[z].amplicon_flux || variants[z].depth_flag || variants[z].outside_freq_range || variants[z].qual_flag || variants[z].del_flag) continue;
        //this pos was flagged as poorly assigned
        if(tmp != assignment_flagged.end() && k == pos_idxs[j]){
          //technically this could use work as it's repetitive
          //std::cerr << variants[z].position << " " << variants[z].freq << " " << assigned[j] << std::endl;
          variants[z].vague_assignment = true; 
          variants[z].cluster_assigned = assigned[j];
          break;
        }
        if(k == pos_idxs[j] && tmp == assignment_flagged.end()){
          uint32_t idx = assigned[j];
          variants[z].cluster_assigned = idx;
        }
        k++;
      }
    }
  }
}
void go(uint32_t offset, uint32_t k, std::vector<double> means, std::vector<double> combination, std::vector<std::vector<double>> &combos, double error) {
  //generates all the combinations
  if (k == 0) {
    double sum_of_elems = std::accumulate(combination.begin(), combination.end(), 0.0);
    //std::cerr << sum_of_elems << std::endl;
    if(sum_of_elems < 1+error && sum_of_elems > 1-error){
      combos.push_back(combination);
    }
    return;
  }
  for (uint32_t i = offset; i <= means.size() - k; ++i) {
    combination.push_back(means[i]);
    go(i+1, k-1, means, combination, combos, error);
    combination.pop_back();
  }
}
void solve_solution_sets(std::vector<double> means, uint32_t n){
  if(means.size() < n){
    return;
  }
  //generate the combinations
  std::vector<double> combination; //placeholder to use in function
  std::vector<std::vector<double>> combos;
  double error = 0.05;
  go(0, n, means, combination, combos, error);
  for(uint32_t i=0; i < combos.size(); i++){
    for(uint32_t j=0; j < combos[i].size(); j++){
      //std::cerr << combos[i][j] << " ";
    }
    //std::cerr << std::endl;
  }

}

void determine_low_prob_positions(std::vector<variant> &variants){
  std::vector<uint32_t> low_prob_pos;
  for(uint32_t i=0; i< variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
      if(variants[i].cluster_assigned == -1) continue;
      double prob = exp(variants[i].probabilities[variants[i].cluster_assigned]);
      if(prob < 0.00001){
        low_prob_pos.push_back(variants[i].position);
      }
    }
  }
  //flag these positions with low confidence in the variants
  for(uint32_t i=0; i< variants.size(); i++){
    std::vector<uint32_t>::iterator it = std::find(low_prob_pos.begin(), low_prob_pos.end(), variants[i].position); 
    if(it != low_prob_pos.end()){
      variants[i].low_prob_flag = true;
    }
  }
}

void determine_outlier_variants(std::vector<variant> &variants, uint32_t n){ 
  std::vector<std::vector<double>> clusters;
  for(uint32_t i=0; i < n; i++){
    std::vector<double> tmp;
    clusters.push_back(tmp);
  }
  for(uint32_t z=0; z < variants.size(); z++){
    if(variants[z].cluster_assigned > -1){
      clusters[variants[z].cluster_assigned].push_back(variants[z].freq);
    }
  }
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  std::vector<double> cluster_means; 
  //look at spread of variants assigned to each cluster to determine outliers
  for(uint32_t i=0; i < clusters.size(); i++){
    double sum = std::accumulate(clusters[i].begin(), clusters[i].end(), 0.0);
    double mean = sum / clusters[i].size();
    double sq_sum = std::inner_product(clusters[i].begin(), clusters[i].end(), clusters[i].begin(), 0.0);
    double stdev = std::sqrt(sq_sum / clusters[i].size() - mean * mean);
    double upper_bound = (stdev*3) + mean;
    double lower_bound = mean - (stdev*3);
    cluster_means.push_back(mean);
    lower_bounds.push_back(lower_bound);
    upper_bounds.push_back(upper_bound);
    //std::cerr << "lower " << lower_bound << " upper " << upper_bound << std::endl;
    //std::cerr << "mean " << mean << " std dev " << stdev << std::endl;
  }
  
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].cluster_assigned > -1){
      double lb = lower_bounds[variants[i].cluster_assigned];
      double ub = upper_bounds[variants[i].cluster_assigned];
      double cm = cluster_means[variants[i].cluster_assigned];
      if((variants[i].freq < lb || variants[i].freq > ub) && std::abs(cm - variants[i].freq) > 0.03){
        variants[i].cluster_outlier = true;
      }
    }
  }
}

void split(const std::string &s, char delim, std::vector<std::string> &elems){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<uint32_t> find_low_quality_positions(std::string filename, uint32_t depth_cutoff, float lower_bound, float upper_bound, float quality_threshold, uint32_t round_val){
  std::vector<uint32_t> low_quality_positions;
  std::ifstream infile(filename);
  std::string line;
  uint32_t count = 0;
  uint32_t pos = 0;
  float qual = 0;
  float depth = 0;
  float freq = 0;
  while (std::getline(infile, line)) {
    if(count == 0) {
      count += 1;
      continue;
    }
    std::vector<std::string> row_values;
    //split the line by delimiter
    split(line, '\t', row_values);
    depth = std::stoi(row_values[2]);
    pos = std::stoi(row_values[0]);
    qual = std::stof(row_values[4]);
    freq = std::stof(row_values[3]);
    float multiplier = pow(10, round_val);
    freq = round(freq * multiplier) / multiplier;
    std::string nuc = row_values[1];
    if(1/freq * depth< depth_cutoff || freq < lower_bound || freq > upper_bound || is_substring(nuc, "+") || is_substring(nuc, "-")) continue;
    if (qual < quality_threshold){ 
      //std::cerr << "position " << pos << " qual " << qual << " depth " << depth <<  " freq " << freq << " nuc " << nuc << std::endl;
      low_quality_positions.push_back(pos);
    }
  }
  return(low_quality_positions);
}

std::vector<uint32_t> find_deletion_positions(std::string filename, uint32_t depth_cutoff, float lower_bound, float upper_bound, uint32_t round_val){

  std::vector<uint32_t> deletion_positions;
  std::ifstream infile(filename);
  std::string line;
  uint32_t count = 0;
  float freq = 0;
  uint32_t depth = 0;
  uint32_t pos = 0;
 
  while (std::getline(infile, line)) {
    if(count == 0) {
      count += 1;
      continue;
    }
    std::vector<std::string> row_values;
    //split the line by delimiter
    split(line, '\t', row_values);
    depth = std::stoi(row_values[2]);
    pos = std::stoi(row_values[0]);
    freq = std::stof(row_values[3]);
    float multiplier = pow(10, round_val);
    freq = round(freq * multiplier) / multiplier;
    std::string nuc = row_values[1];
    if(1/freq * depth < depth_cutoff || freq < lower_bound || freq > upper_bound) continue;
    if (is_substring(nuc, "+") || is_substring(nuc, "-")) {
      deletion_positions.push_back(pos);
      //std::cerr << "del pos depth " << depth << " freq " << freq << std::endl; 
    }
  }
  return(deletion_positions);
}


void parse_internal_variants(std::string filename, std::vector<variant> &variants, uint32_t depth_cutoff, float lower_bound, float upper_bound, std::vector<uint32_t> deletion_positions, std::vector<uint32_t> low_quality_positions, uint32_t round_val){
  /*
   * Parses the variants file produced internally by reading bam file line-by-line.
   */
  std::ifstream infile(filename);
  std::string line;
  uint32_t count = 0;
  while (std::getline(infile, line)) {
    if(count == 0) {
      count += 1;
      continue;
    }
    std::vector<std::string> row_values;
    //split the line by delimiter
    uint32_t pos = 0;
    uint32_t depth = 0;
    float freq = 0;
    float qual = 0;
    std::string flag = "";
    std::string is_ref = "";
    split(line, '\t', row_values);

    pos = std::stoi(row_values[0]);
    bool del_pos = std::find(deletion_positions.begin(), deletion_positions.end(), pos) != deletion_positions.end();
    bool lq_pos = std::find(low_quality_positions.begin(), low_quality_positions.end(), pos) != low_quality_positions.end();

    depth = std::stoi(row_values[2]);
    freq = std::stof(row_values[3]);

    float multiplier = pow(10, round_val);
    freq = round(freq * multiplier) / multiplier;
    qual = std::stof(row_values[4]);
    flag = row_values[5];
    is_ref = row_values[6];
   
    variant tmp;
    tmp.position = pos;
    tmp.nuc = row_values[1];
  
    tmp.depth = depth;
    tmp.freq = freq;
    tmp.qual = qual;

    if (flag == "TRUE"){
      tmp.amplicon_flux = true;
    } else {
      tmp.amplicon_flux = false;
    }
    if(1/freq * depth < depth_cutoff){
      tmp.depth_flag = true;
    } else {
      tmp.depth_flag = false;
    } 
    if (freq < lower_bound || freq > upper_bound){
      tmp.outside_freq_range = true;
    } else {
      tmp.outside_freq_range = false;
    }
    if (is_ref == "TRUE"){
      tmp.is_ref = true;
    } else {
      tmp.is_ref = false;
    }
    if(lq_pos){
      tmp.qual_flag = true;
    } else {
      tmp.qual_flag = false;
    }
    if(del_pos){
      tmp.del_flag = true;
    } else {
      tmp.del_flag = false;
    }
    variants.push_back(tmp); 
    count += 1;
  } 
}

int gmm_model(std::string prefix, std::vector<uint32_t> populations_iterate, std::string output_prefix){
  int retval = 0;
  float lower_bound = 0.03;
  float upper_bound = 0.99;
  uint32_t depth_cutoff = 50;
  float quality_threshold = 20;
  uint32_t round_val = 3;
  std::vector<variant> variants;
  std::vector<uint32_t> deletion_positions = find_deletion_positions(prefix, depth_cutoff, lower_bound, upper_bound, round_val);
  std::vector<uint32_t> low_quality_positions = find_low_quality_positions(prefix, depth_cutoff, lower_bound, upper_bound, quality_threshold, round_val);
  
  parse_internal_variants(prefix, variants, depth_cutoff, lower_bound, upper_bound, deletion_positions, low_quality_positions, round_val);
  std::string filename = prefix + ".txt";
  //this whole things needs to be reconfigured
  uint32_t useful_var=0;
  for(uint32_t i=0; i < variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
      useful_var += 1;
    }
  }
  //initialize armadillo dataset and populate with frequency data
  arma::mat data(1, useful_var, arma::fill::zeros);
  //(rows, cols) where each columns is a sample  
  uint32_t count=0;
  for(uint32_t i = 0; i < variants.size(); i++){
    //check if variant should be filtered for first pass model
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
      double tmp = static_cast<double>(variants[i].freq);
      data.col(count) = tmp;
      count += 1;
    }
  }
  std::vector<std::vector<double>> solutions; //straight from the model
  std::vector<std::vector<double>> solutions_recalculated; //solutions without outliers
  std::vector<double> all_aic; //aic for each population
  std::vector<uint32_t> unassigned; //count number of variants can't be assigned
  std::vector<std::vector<double>> stddevs; //save the standard deviations of each cluster
  for(auto n : populations_iterate){
    variants.clear();
    parse_internal_variants(prefix, variants, depth_cutoff, lower_bound, upper_bound, deletion_positions, low_quality_positions, round_val);

    if(n > (useful_var/10)) continue; //this is because it's recommended to have 10 points per gaussian
    //model learning
    arma::gmm_diag model;
    bool status = model.learn(data, n, arma::eucl_dist, arma::random_spread, 10, 15, 1e-10, false);
    if(status == false){
      std::cerr << "gmm model failed" << std::endl;
      continue;
    }
    //get the means of the gaussians
    std::vector<double> means;
    //means with the outliers removed
    //model.means.print("means:");
    //get the probability of each frequency being assigned to each gaussian
    std::vector<std::vector<double>> prob_matrix;
    for(uint32_t i=0; i < n; i++){
      means.push_back((double)model.means[i]);
      arma::rowvec set_likelihood = model.log_p(data.cols(0,useful_var-1), i);
      std::vector<double> tmp;
      for(uint32_t j=0; j < useful_var; j++){
        tmp.push_back((double)set_likelihood[j]);
      }
      prob_matrix.push_back(tmp);
    }
    std::vector<std::vector<double>> tv = transpose_vector(prob_matrix);
    uint32_t j = 0;
    for(uint32_t i=0; i < variants.size(); i++){
      if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
        variants[i].probabilities = tv[j];
        j++;
      }
    }
    //assign variants out based on probability, not taking into account condition of all variants for a pos ~= 1 
    assign_variants_simple(variants, prob_matrix);
  
    //calculate cluster outliers
    determine_outlier_variants(variants, n);

    uint32_t unassigned_var = 0;
    uint32_t outlier_var =0;
    uint32_t low_prob_var = 0;
    double prob_sum = 0;
    std::vector<std::vector<double>> clusters(n); 

    determine_low_prob_positions(variants);

    bool useful = true;
    std::cerr << "useful variants " << useful_var << std::endl;
    for(uint32_t i=0; i < variants.size(); i++){
      if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range && !variants[i].qual_flag && !variants[i].del_flag){
        if(variants[i].cluster_assigned == -1){
          unassigned_var += 1;
         }
        //log likelihood for each point
        double prob = variants[i].probabilities[variants[i].cluster_assigned];
        prob_sum += prob;
        if(variants[i].cluster_assigned != -1){
          double var = (double)variants[i].freq;
          clusters[variants[i].cluster_assigned].push_back(var);
        }
      }
    }
    std::vector<double> recalculated_means;
    std::vector<double> tmp_stddev;
    useful = true;
    for(uint32_t i=0; i < clusters.size(); i++){
      if(clusters[i].size() == 0){
        useful = false;
        break;
      }
      std::cerr << "i " << i << std::endl;
      double sum = std::accumulate(clusters[i].begin(), clusters[i].end(), 0.0);
      double mean = sum / clusters[i].size();
      double sq_sum = std::inner_product(clusters[i].begin(), clusters[i].end(), clusters[i].begin(), 0.0);
      double stdev = std::sqrt(sq_sum / clusters[i].size() - mean * mean);
      //recalculated_means.push_back(mean);  
      tmp_stddev.push_back(stdev);
      /*for(auto x : clusters[i]){
        std::cerr << x << " ";
      }
      std::cerr << std::endl;*/
      recalculated_means.push_back(mean);
      std::cerr << "mean " << mean << " std dev " << stdev << std::endl;
      //std::cerr << "median " << calculate_median(clusters[i]) << std::endl;
    }
    if(!useful) continue;
    solutions.push_back(means);
    stddevs.push_back(tmp_stddev);
    solutions_recalculated.push_back(recalculated_means);
    unassigned.push_back(unassigned_var+outlier_var+low_prob_var);
    //std::cerr << useful_var << " probability " << prob_sum <<  " avg prob " << prob_sum / useful_var << std::endl;
    double aic = (2 * (double)n) - (2 * prob_sum / useful_var);
    all_aic.push_back(aic);
    std::cerr << "aic " << aic << std::endl;
    std::cerr << "unassigned variants " << unassigned_var << " outlier var " << outlier_var << " low prob var " << low_prob_var << std::endl;
    //std::cerr << "\n";
    //draw the actual threshold
    //double threshold = calculate_cluster_bounds(variants, n);
    //std::cerr << "n " << n << " threshold " << threshold << std::endl;
    //call the consensus sequence
    //model.save("my_model.gmm");
  }
  auto min_idx_iter = std::min_element(all_aic.begin(), all_aic.end()); 
  int index = std::distance(all_aic.begin(), min_idx_iter);
  for(uint32_t i=1; i < all_aic.size(); i++){
    if(all_aic[i-1] < all_aic[i]){
      //std::cerr << i-1 << std::endl;
      index = (int)i-1;
      break;
    }
  }
  //here we now define the two criteria in which we eliminate things as being contaminated
  //1. the solution is not solveable for identity reasons
  //2. the clustered solution had too many outliers
  //3. unexplained peaks
  std::vector<double> means = solutions_recalculated[index];
  std::vector<double> stddev = stddevs[index];
  uint32_t unassigned_var = unassigned[index];

  //here we make sure that we have no unexplained peaks  
  std::vector<std::vector<double>> viable_solutions = solve_possible_solutions(means); 
  if(viable_solutions.size() == 0){
    std::cerr << "no viable solutions" << std::endl;
    retval = -1;
  }

  std::vector<std::vector<double>> kept_solutions;
  if(viable_solutions.size() > 0){  
    kept_solutions = deduplicate_solutions(viable_solutions);
    if(kept_solutions.size() == 0){
      std::cerr << "no solutions after deduplication" << std::endl;
      retval = -1;
    }
  } 
  //check to see if this is solveable
  std::vector<uint32_t> sizes;
  for(auto &vec : kept_solutions){
    std::sort(vec.begin(), vec.end(), std::greater<double>());
    sizes.push_back(vec.size());
    /*for(auto x : vec){
      std::cerr << x << " ";
    }
    std::cerr << std::endl;*/
  } 
  uint32_t max_solution=0;
  if(sizes.size() > 0){
    //here we need each solution to be a set of the other for it to be solveable
    max_solution = *std::max_element(sizes.begin(), sizes.end());  
    max_solution = 1;
  }
  bool solveable = true;
  for(uint32_t i=0; i < max_solution; i++){
    std::vector<double> values;
    for(auto vec : kept_solutions){
      if(i < vec.size()){
        values.push_back(vec[i]);
      }
    }
    //bool solveable = true;
    double first_val = values[0];
    for(auto x : values){
      if(x != first_val){ 
        solveable = false;
        break;
      }
    }
  }
  std::string stddev_string = "[";
  std::string peaks_string = "[";
  for(uint32_t i=0; i < means.size(); i++){
    if(i != 0){
      stddev_string += ",";
      peaks_string += ",";
    }
    stddev_string += std::to_string(stddev[i]);
    peaks_string += std::to_string(means[i]);
  }
  stddev_string += "]";
  peaks_string += "]";
  std::string means_string = "[";
  for(uint32_t j=0; j < kept_solutions.size(); j++){
    means_string += "[";
    for(uint32_t i=0; i < kept_solutions[j].size(); i++){
      if(i != 0) means_string += ",";
      means_string += std::to_string(kept_solutions[j][i]);
    }
    if(j != kept_solutions.size()-1){
      means_string += "],";
    } else {
      means_string += "]";
    }
  }
  means_string += "]";
  std::ofstream file;  
  file.open(output_prefix + ".txt", std::ios::trunc);
  file << "n\tmeans\tstdev\tsolutions\taic\tnum_solutions\tidentity_issue\tunassigned_var\tuseful_var\n";
  file << std::to_string(populations_iterate[index]) << "\t";
  file << peaks_string << "\t";
  file << stddev_string << "\t";
  file << means_string << "\t";
  file << std::to_string(all_aic[index]) << "\t";
  file << std::to_string(kept_solutions.size()) << "\t";
  if(solveable){
    file << "false\t";
  } else {
    file << "true\t";
  }
  file << std::to_string(unassigned_var) << "\t";
  file << std::to_string(useful_var) << "\n";
  file.close();

  cluster_consensus();

  return(retval);
}
