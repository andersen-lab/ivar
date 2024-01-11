#include "./include/armadillo"
#include <fstream>
#include <cmath>

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

std::vector<std::vector<double>> assign_variants_simple(std::vector<float> filtered_frequencies, std::vector<uint32_t> filtered_positions, std::vector<std::vector<double>> prob_matrix){
  /*
   * Here we assign variants to the highest probability group on a per position basis making sure two variants at the same poistion don't end up in the same group.
   */
  uint32_t n = prob_matrix.size();
  std::vector<std::vector<double>> clusters;
  //populate empty cluster vectors
  for(uint32_t i=0; i < n; i++){
    std::vector<double> tmp;
    clusters.push_back(tmp);
  }
  //find the unique positions
  std::vector<uint32_t> unique_pos;
  for(uint32_t i = 0; i < filtered_positions.size(); i++){
    if (std::find(unique_pos.begin(), unique_pos.end(), filtered_positions[i]) == unique_pos.end()) {
      unique_pos.push_back(filtered_positions[i]);
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
    //all locations in the prob matrix for this position
    for(uint32_t j = 0; j < filtered_positions.size(); j++){
      if(filtered_positions[j] == unique_pos[i]){
        pos_idxs.push_back(j);
        std::vector<double> tmp;
        for(uint32_t k=0; k < n; k++){
          tmp.push_back(prob_matrix[k][j]);
        }
        tmp_prob.push_back(tmp);
      }
    }
    //assign variants based on most probable position-wise
    std::vector<uint32_t> assigned = calculate_joint_probabilities(tmp_prob, possible_permutations);
    //make sure the assignment is concrete
    std::vector<uint32_t> assignment_flagged = compare_cluster_assignment(tmp_prob, assigned);
    for(uint32_t j=0; j < pos_idxs.size(); j++){
      if(std::find(assignment_flagged.begin(), assignment_flagged.end(), j) != assignment_flagged.end()) { 
        continue;
      }
      uint32_t idx = assigned[j];
      clusters[idx].push_back(filtered_frequencies[pos_idxs[j]]);
    }
  }
  return(clusters);
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

void determine_outlier_variants(std::vector<std::vector<double>> clusters, std::vector<uint32_t> &outlier_idxs, std::vector<double> means){

  //look at spread of variants assigned to each cluster to determine outliers
  /*for(uint32_t i=0; i < clusters.size(); i++){
    for(uint32_t j=0; j < clusters[i].size(); j++){
      std::cerr << clusters[i][j] << " ";
    }
    std::cerr << std::endl;
    auto const Q1 = clusters[i].size() / 4;
    auto const Q2 = clusters[i].size() / 2;
    auto const Q3 = Q1 + Q2;

    std::nth_element(clusters[i].begin(), clusters[i].begin() + Q1, clusters[i].end());
    std::nth_element(clusters[i].begin() + Q1 + 1, clusters[i].begin() + Q2, clusters[i].end());
    std::nth_element(clusters[i].begin() + Q2 + 1, clusters[i].begin() + Q3, clusters[i].end());
    std::cerr << Q1 << " " << clusters[i][Q1] << std::endl;
    std::cerr << Q3 << " " << clusters[i][Q3] << std::endl;
    double sum = std::accumulate(clusters[i].begin(), clusters[i].end(), 0.0);
    double mean = sum / clusters[i].size();
    double sq_sum = std::inner_product(clusters[i].begin(), clusters[i].end(), clusters[i].begin(), 0.0);
    double stdev = std::sqrt(sq_sum / clusters[i].size() - mean * mean);

    double upper_bound = (stdev*3) + mean;
    double lower_bound = mean - (stdev*3);
    std::cerr << "lower " << lower_bound << " upper " << upper_bound << std::endl;
    std::cerr << "mean " << mean << " std dev " << stdev << std::endl;
  }*/

}

void split(const std::string &s, char delim, std::vector<std::string> &elems){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

void parse_internal_variants(std::string filename, std::vector<float> &frequencies, std::vector<uint32_t> &positions, std::vector<uint32_t> &depths, std::vector<std::string> &flagged, std::vector<std::string> &nucs){
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
    split(line, '\t', row_values);
    pos = std::stoi(row_values[0]);
    nucs.push_back(row_values[1]);
    depth = std::stoi(row_values[2]);
    freq = std::stof(row_values[3]);
    flag = row_values[5];
    frequencies.push_back(freq);
    positions.push_back(pos);
    depths.push_back(depth);
    flagged.push_back(flag);
    count += 1;
  } 
}

int gmm_model(std::string prefix){
  arma::gmm_diag model;
  int retval = 0;
  float lower_bound = 0.03;
  float upper_bound = 0.97;

  std::string filename = prefix + ".txt";

  std::vector<float> frequencies;
  std::vector<uint32_t> positions;
  std::vector<uint32_t> depths;
  std::vector<std::string> nucs;
  std::vector<std::string> flagged;
  parse_internal_variants(filename, frequencies, positions, depths, flagged, nucs);

  //the n min and n max represent the range of n values
  uint32_t n = 3;

  //filter out the data we won't use in our initial model
  std::vector<float> filtered_frequencies;
  std::vector<uint32_t> filtered_positions;
  for(uint32_t i=0; i < frequencies.size(); i++){
    if(depths[i] > 10 && flagged[i] == "FALSE" && frequencies[i] > lower_bound && frequencies[i] < upper_bound){
      filtered_frequencies.push_back(frequencies[i]);  
      filtered_positions.push_back(positions[i]);
    }
  }

  //initialize armadillo dataset and populate with frequency data
  //(rows, cols) where each columns is a sample
  arma::mat data(1, filtered_frequencies.size(), arma::fill::zeros); 
  for(uint32_t i = 0; i < filtered_frequencies.size(); i++){
    double tmp = static_cast<double>(filtered_frequencies[i]);
    data.col(i) = tmp;
  }
  //model learning
  bool status = model.learn(data, n, arma::eucl_dist, arma::random_subset, 10, 5, 1e-10, true);
  if(status == false){
    std::cerr << "gmm model failed" << std::endl;
  }
  //get the means of the gaussians
  std::vector<double> means;
  
  //TEST LINES
  for(uint32_t i=0; i < filtered_frequencies.size(); i++){
    std::cerr << filtered_frequencies[i] << " " << filtered_positions[i] << " " << model.log_p(data.col(i), 0) << " " << model.log_p(data.col(i), 1) << " " << model.log_p(data.col(i), 2) << std::endl;
  }

  //get the probability of each frequency being assigned to each gaussian
  std::vector<std::vector<double>> prob_matrix;
  for(uint32_t i=0; i < n; i++){
    //means.push_back((double)model.means[i]);
    arma::rowvec set_likelihood = model.log_p(data.cols(0,filtered_frequencies.size()-1), i);
    std::vector<double> tmp;
    for(uint32_t j=0; j < filtered_frequencies.size(); j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }
    
  //assign variants out based on probability, not taking into account condition of all variants for a pos ~= 1 
  std::vector<std::vector<double>> clusters = assign_variants_simple(filtered_frequencies, filtered_positions, prob_matrix);

  //calculate cluster means
  for(uint32_t i=0; i < clusters.size(); i++){
    double m = accumulate(clusters[i].begin(), clusters[i].end(), 0.0) / clusters[i].size();
    std::cerr << m << std::endl;
    means.push_back(m);
  }
  //based on the probabilities, determine which points cannot be concretely assigned
  std::vector<uint32_t> outlier_idxs; 
  determine_outlier_variants(clusters, outlier_idxs, means);
  exit(1);

  //solve for solutions adding to ~1
  //solve_solution_sets(means, n);
  
  //places where we'll call an N value
  //std::vector<uint32_t> unassigned_points;
  //model.save("my_model.gmm");
  return(retval);
}
