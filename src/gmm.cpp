#include "./include/armadillo"
#include "gmm.h"
#include <fstream>
#include <cmath>

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

std::vector<uint32_t> estimate_populations(std::vector<variant> variants){
  //TODO generalize
  uint32_t num_var = 0;
  for(uint32_t i=0; i < variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && variants[i].freq > 0.01 && !variants[i].is_ref){
      num_var += 1;
    }
  }
  //this isn't exactly what we do for ww so needs to be tested thoroughly
  float evolutionary_rate = 0.001;
  float time_elapsed = 0.83; //in years
  float ref = 29003;
  float beta = 1;
  float percent_possible_mutations = (num_var * beta) / (ref * 3);
  float expected_nt_sub_per_site = evolutionary_rate * time_elapsed;
  float num_populations = percent_possible_mutations / expected_nt_sub_per_site;

  std::vector<uint32_t> population_bounds;
  //can't have 0 or fewer populations present
  uint32_t n = round(num_populations);
  if (n-1 <= 0) {
    population_bounds.push_back(1);
  } else {
    population_bounds.push_back(n-1);
  }
  population_bounds.push_back(n+1);
  return(population_bounds);
}

uint32_t count_useful_variants(std::vector<variant> variants){
  uint32_t count = 0;
  //determine the number of variants useful for modeling
  for(uint32_t i=0; i< variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range){
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

std::vector<uint32_t>  compare_cluster_assignment(std::vector<std::vector<double>> prob_matrix, std::vector<uint32_t> assigned){
  double threshold = 3;
  std::vector<uint32_t> flagged_idx;
  for(uint32_t i=0; i < prob_matrix.size(); i++){
    double assigned_prob = prob_matrix[i][assigned[i]];
    std::vector<double> tmp = prob_matrix[i];
    std::sort(tmp.begin(), tmp.end(), std::greater<double>());
    for(uint32_t j=0; j < tmp.size(); j++){
      if(tmp[j] >= assigned_prob) continue;
        if(exp(tmp[j]) * threshold > exp(assigned_prob)){
          flagged_idx.push_back(i);
        }
        break;
    }
  }
  return(flagged_idx);
}

std::vector<uint32_t>  calculate_joint_probabilities(std::vector<std::vector<double>> prob_matrix, std::vector<std::vector<uint32_t>> permutations){
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
        score = -1000;
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
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range){
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
      if(variants[k].amplicon_flux || variants[k].depth_flag || variants[k].outside_freq_range) continue;
      
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
        if(variants[z].amplicon_flux || variants[z].depth_flag || variants[z].outside_freq_range) continue;
        //this pos was flagged as poorly assigned
        if(tmp != assignment_flagged.end()){
          //technically this could use work as it's repetitive
          //std::cerr << variants[z].position << " " << variants[z].freq << std::endl;
          variants[z].vague_assignment = true;
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
    std::cerr << sum_of_elems << std::endl;
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
      std::cerr << combos[i][j] << " ";
    }
    std::cerr << std::endl;
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
  //look at spread of variants assigned to each cluster to determine outliers
  for(uint32_t i=0; i < clusters.size(); i++){
    double sum = std::accumulate(clusters[i].begin(), clusters[i].end(), 0.0);
    double mean = sum / clusters[i].size();
    double sq_sum = std::inner_product(clusters[i].begin(), clusters[i].end(), clusters[i].begin(), 0.0);
    double stdev = std::sqrt(sq_sum / clusters[i].size() - mean * mean);
    double upper_bound = (stdev*3) + mean;
    double lower_bound = mean - (stdev*3);
    lower_bounds.push_back(lower_bound);
    upper_bounds.push_back(upper_bound);
    //std::cerr << "lower " << lower_bound << " upper " << upper_bound << std::endl;
    //std::cerr << "mean " << mean << " std dev " << stdev << std::endl;
  }
  
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].cluster_assigned > -1){
      double lb = lower_bounds[variants[i].cluster_assigned];
      double ub = upper_bounds[variants[i].cluster_assigned];
      if(variants[i].freq < lb || variants[i].freq > ub){
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

void parse_internal_variants(std::string filename, std::vector<variant> &variants, uint32_t depth_cutoff, float lower_bound, float upper_bound){
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
    std::string flag = "";
    std::string is_ref = "";
    split(line, '\t', row_values);


    pos = std::stoi(row_values[0]);
    depth = std::stoi(row_values[2]);
    freq = std::stof(row_values[3]);
    flag = row_values[5];
    is_ref = row_values[6];
   
    variant tmp;
    tmp.position = pos;
    tmp.nuc = row_values[1];
    tmp.depth = depth;
    tmp.freq = freq;
    if (flag == "TRUE"){
      tmp.amplicon_flux = true;
    } else {
      tmp.amplicon_flux = false;
    }
    if (depth < depth_cutoff){
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
    variants.push_back(tmp); 
    count += 1;
  } 
}

int gmm_model(std::string prefix){
  arma::gmm_diag model;
  int retval = 0;
  float lower_bound = 0.03;
  float upper_bound = 0.97;
  uint32_t depth_cutoff = 10;


  std::string filename = prefix + ".txt";
  std::vector<variant> variants;

  parse_internal_variants(filename, variants, depth_cutoff, lower_bound, upper_bound);
  uint32_t useful_var = count_useful_variants(variants); 
  std::vector<uint32_t> population_bounds = estimate_populations(variants);

  //the n min and n max represent the range of n values
  uint32_t n_min = population_bounds[0];
  uint32_t n_max = population_bounds[1];
  uint32_t n = 2; 
  //the number of variants we will acutally use when modeling
  std::cerr << "n min " << n_min << " n max " << n_max << std::endl;
 
  //sample is too complex and we flag as a contaminant
  //TODO make this bit better
  if(n_max > 3){
    retval = 1;
    return(retval);
  }

  //initialize armadillo dataset and populate with frequency data
  //(rows, cols) where each columns is a sample
  
  arma::mat data(1, useful_var, arma::fill::zeros);
  std::cerr << "useful var " << useful_var << std::endl; 
  uint32_t count = 0;
  for(uint32_t i = 0; i < variants.size(); i++){
    //check if variant should be filtered for first pass model
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range){
      double tmp = static_cast<double>(variants[i].freq);
      data.col(count) = tmp;
      count += 1;
    }
  }
  
  //model learning
  bool status = model.learn(data, n, arma::eucl_dist, arma::random_subset, 10, 5, 1e-10, true);
  if(status == false){
    std::cerr << "gmm model failed" << std::endl;
  }
  //get the means of the gaussians
  std::vector<double> means;
  model.means.print("means:");
  uint32_t j=0;
  for(uint32_t i=0; i < variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && !variants[i].outside_freq_range){
      //std::cerr << variants[i].freq << " " << model.log_p(data.col(j), 0) << " " << model.log_p(data.col(j), 1) << " " << model.log_p(data.col(j), 2) << std::endl;
      std::cerr << variants[i].freq << " " << model.log_p(data.col(j), 0) << " " << model.log_p(data.col(j), 1) << std::endl;
      j++;
    }
  }  
  
  //get the probability of each frequency being assigned to each gaussian
  std::vector<std::vector<double>> prob_matrix;
  for(uint32_t i=0; i < n; i++){
    //means.push_back((double)model.means[i]);
    arma::rowvec set_likelihood = model.log_p(data.cols(0,useful_var-1), i);
    std::vector<double> tmp;
    for(uint32_t j=0; j < useful_var; j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }
   
  //assign variants out based on probability, not taking into account condition of all variants for a pos ~= 1 
  assign_variants_simple(variants, prob_matrix);
   
  //calculate cluster outliers
  determine_outlier_variants(variants, n);

  //draw the actual threshold
  double threshold = calculate_cluster_bounds(variants, n);
  std::cerr << "threshold " << threshold << std::endl;
  //call the consensus sequence

  //model.save("my_model.gmm");

  return(retval);
}
