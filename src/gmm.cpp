#include "./include/armadillo"
#include "gmm.h"
#include "call_consensus_clustering.h"
#include "estimate_error.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>

double calculate_range(const std::vector<double>& vec) {
    if (vec.empty()) {
        throw std::invalid_argument("Vector must not be empty.");
    }

    // Find the maximum and minimum values in the vector
    double max_value = *std::max_element(vec.begin(), vec.end());
    double min_value = *std::min_element(vec.begin(), vec.end());

    // Calculate the range
    return max_value - min_value;
}

double calculate_variance(std::vector<double> arr, uint32_t size) {
    double mean = 0.0;
    double variance = 0.0;

    // Calculate the mean
    for (uint32_t i = 0; i < size; ++i) {
        mean += arr[i];
    }
    mean /= size;

    // Calculate the variance
    for (uint32_t i = 0; i < size; ++i) {
        variance += std::pow(arr[i] - mean, 2);
    }
    variance /= size; // For population variance
    // variance /= (size - 1); // Uncomment this line for sample variance

    return variance;
}

// Function to calculate the mean of a vector
double calculate_mean(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0f;
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0f);
    return sum / data.size();
}

double calculate_median(std::vector<double> &data) {
    std::sort(data.begin(), data.end());
    size_t size = data.size();

    if (size % 2 == 0) {
        return (data[size / 2 - 1] + data[size / 2]) / 2.0;
    } else {
        return data[size / 2];
    }
}

double calculate_stddev(const std::vector<double>& data, double mean) {
    double sum = 0;
    for (double val : data) {
        sum += std::pow(val - mean, 2);
    }
    return std::sqrt(sum / data.size());
}

int count_outliers_z_score(const std::vector<double>& data) {
    double mean = calculate_mean(data);
    double stddev = calculate_stddev(data, mean);

    int outlier_count = 0;
    for (double value : data) {
        double z_score = (value - mean) / stddev;
        if (std::abs(z_score) > 3) { // Typically consider |Z| > 3 as an outlier
            outlier_count++;
        }
    }

    return outlier_count;
}

void calculate_gapped_frequency(std::vector<variant> &variants, double universal_cluster, double noise_cluster){
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].pos_del_flag){
      double total_depth = 0;
      for(uint32_t j=0; j < variants.size(); j++){
        if(variants[i].position == variants[j].position){
          total_depth += variants[j].depth;
        }
      }
      variants[i].gapped_freq = variants[i].depth / total_depth;
      variants[i].transformed_gap_freq = log(variants[i].gapped_freq / (1 - variants[i].gapped_freq));
      if(variants[i].gapped_freq > universal_cluster || variants[i].gapped_freq < noise_cluster){
        variants[i].gap_outside_freq_range = true;
      }
    }
  }
}

// Function to calculate the interquartile range (IQR)
double calculate_interquartile_range(std::vector<double> data) {
    std::sort(data.begin(), data.end()); // Sort the data
    
    size_t n = data.size();
    size_t mid = n / 2;

    // Split data into lower half (Q1 range) and upper half (Q3 range)
    std::vector<double> lower_half(data.begin(), data.begin() + mid);
    std::vector<double> upper_half(data.begin() + (n % 2 == 0 ? mid : mid + 1), data.end());

    // Calculate Q1 (median of lower half) and Q3 (median of upper half)
    double q1 = calculate_median(lower_half);
    double q3 = calculate_median(upper_half);

    // Return the interquartile range
    return q3 - q1;
}


uint32_t find_max_frequency_count(const std::vector<uint32_t>& nums) {
    std::unordered_map<uint32_t, uint32_t> frequency_map;
    uint32_t max_count = 0;

    // Count frequency of each number
    for (const auto& num : nums) {
        frequency_map[num]++;
        if (frequency_map[num] > max_count) {
            max_count = frequency_map[num];
        }
    }

    return max_count;
}

double calculate_mad(const std::vector<double>& data, double mean){
    //Calculate the sum of absolute deviations
    double absDevSum = 0.0;
    for (double value : data) {
        absDevSum += std::abs(value - mean);
    }

    // Calculate and return the MAD
    return absDevSum / data.size();
}

void assign_clusters(std::vector<variant> &variants, gaussian_mixture_model gmodel, uint32_t lower_n){
  std::vector<std::vector<double>> tv = transpose_vector(gmodel.prob_matrix);
  uint32_t j = 0;
  for(uint32_t i=0; i < variants.size(); i++){
    variants[i].probabilities = tv[j];
    j++;
  }
  uint32_t index = smallest_value_index(gmodel.means);
  assign_variants_simple(variants, gmodel.prob_matrix, index, gmodel.means, lower_n);

}

uint32_t count_repeated_values(const std::vector<double>& vec) {
    std::unordered_map<double, int> frequency_map;
    for (double value : vec) {
        frequency_map[value]++;
    }

    std::vector<double> repeats;
    uint32_t count=0;
    for (const auto& pair : frequency_map) {
      if (pair.second > 1) {
        repeats.push_back(pair.first);
        count += pair.second;
        count -= 1;
      }
    }
    return(count);
}

//function used for production
gaussian_mixture_model retrain_model(uint32_t n, arma::mat data, std::vector<variant> variants, uint32_t lower_n, double var_floor){

  gaussian_mixture_model gmodel; 
  gmodel.n = n;
  // Matrix to store the centroids
  arma::mat centroids;
  arma::mat initial_means(1, n, arma::fill::zeros);

  arma::mat cov (1, n, arma::fill::zeros);
  std::vector<double> total_distances;
  std::vector<std::vector<double>> all_centroids;
 
  for(uint32_t j=0; j < 5; j++){
    //std::cerr << "iteration j " << j << std::endl;
    bool status2 = arma::kmeans(centroids, data, n, arma::random_subset, 15, false);
    if(!status2) continue;

    double total_dist = 0;
    std::vector<std::vector<double>> clusters(n);
    for(auto point : data){
      //using std::min_element to find the closest element
      auto closest_it = std::min_element(centroids.begin(), centroids.end(),
        [point](double a, double b) {
            return std::abs(a - point) < std::abs(b - point);
        });

       uint32_t index = std::distance(centroids.begin(), closest_it);
       //std::cerr << point << " " << centroids[index] << std::endl;
       //clusters[index].push_back(point);
       total_dist += std::abs(point-centroids[index]);
    }    
    /*std::vector<double> vars;
    for(auto c : clusters){
      double var = calculate_variance(c, c.size());
      vars.push_back(var);
    }
    double range = calculate_range(vars);
    total_distances.push_back(range);
    std::cerr << range << std::endl;*/
    std::vector<double> tmp;
    for(auto c : centroids){
      std::cerr << c << " ";
      tmp.push_back((double)c);
    }
    std::cerr << "\n";
    all_centroids.push_back(tmp);
    total_distances.push_back(total_dist);
  }
  
  uint32_t i=0;
  auto min_it = std::min_element(total_distances.begin(), total_distances.end());
  uint32_t index = std::distance(total_distances.begin(), min_it);
  std::vector<double> centroid_vec = all_centroids[index];
  
  for(uint32_t c=0; c < centroid_vec.size(); c++){
    //std::cerr << 0.5 + (1.0 / 20) * std::log(centroid_vec[c] / (1.0 - centroid_vec[c])) << std::endl;
    std::cerr << "initial means " << centroid_vec[c] << std::endl;
    initial_means.col(i) = (double)centroid_vec[c];
    //initial_means.col(i) = 0.5 + (1.0 / 1) * std::log(centroid_vec[c] / (1.0 - centroid_vec[c]));
    cov.col(i) = 0.005;
    ++i;
  }
  //original had this a 0.001, then 0.8, then 0.01
  arma::gmm_diag model;
  model.reset(1, n);
  model.set_means(initial_means);
  model.set_dcovs(cov);
  bool status = model.learn(data, n, arma::eucl_dist, arma::keep_existing, 1, 10, var_floor, false);
  if(!status){
    std::cerr << "model failed to converge" << std::endl;
  }
  
  for(auto m : model.means){
    std::cerr << "retrain means " << m << std::endl;
  }
  std::vector<double> means;
  std::vector<double> hefts;
  std::vector<double> dcovs;
  arma::mat mean_fill2 (1, n, arma::fill::zeros);

  for(uint32_t i=0; i < model.means.size(); i++){
    double m = (double)model.means[i];
    double factor = std::pow(10.0, 2);
    double rounded = std::round(m * factor) / factor;
    mean_fill2.col(i) = rounded;
    means.push_back(rounded);
  }
  model.set_means(mean_fill2);
 
  std::vector<std::vector<double>> prob_matrix;
  std::vector<double> tmp;
  for(uint32_t i=0; i < n; i++){
    arma::rowvec set_likelihood = model.log_p(data, i);
    tmp.clear();
    for(uint32_t j=0; j < data.n_cols; j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }
  gmodel.dcovs = dcovs;
  gmodel.prob_matrix = prob_matrix;
  gmodel.means = means;
  gmodel.hefts = hefts;
  gmodel.model = model;

  return(gmodel);
}

gaussian_mixture_model retrain_model_seeded(uint32_t n, arma::mat data, std::vector<variant> variants, uint32_t lower_n, std::vector<double> centroids){
  gaussian_mixture_model gmodel; 
  gmodel.n = n;

  arma::mat initial_means(1, n, arma::fill::zeros);
  arma::mat cov (1, n, arma::fill::zeros);

  for(uint32_t i=0; i < centroids.size(); i++){
    initial_means.col(i) = centroids[i];
    cov.col(i) = 0.001;
  } 

  double var_floor = 0.001;
  arma::gmm_diag model;
  model.reset(1, n);
  model.set_means(initial_means);
  model.set_dcovs(cov);

  bool status = model.learn(data, n, arma::eucl_dist, arma::keep_existing, 1, 10, var_floor, false);
  if(!status){
    std::cerr << "model failed to converge" << std::endl;
  }

  std::vector<std::vector<double>> prob_matrix;
  std::vector<double> tmp;
  for(uint32_t i=0; i < n; i++){
    arma::rowvec set_likelihood = model.log_p(data, i);
    tmp.clear();
    for(uint32_t j=0; j < data.n_cols; j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }

  std::vector<double> means;
  for(uint32_t i=0; i < model.means.size(); i++){
    double m = (double)model.means[i];
    std::cerr << m << std::endl;
    double factor = std::pow(10.0, 2);
    double rounded = std::round(m * factor) / factor;
    means.push_back(rounded);
  }
  gmodel.prob_matrix = prob_matrix;
  gmodel.means = means;
  gmodel.model = model;

  return(gmodel);
}

double calculate_distance(double point, double mean) {
  return std::abs(point - mean);
}

int find_closest_mean_index(double data_point, const std::vector<double>& means) {
    // Find the index of the closest mean to the data point
    int closest_index = 0;
    double min_distance = std::numeric_limits<double>::infinity();

    for (uint32_t i = 0; i < means.size(); ++i) {
        double distance = calculate_distance(data_point, means[i]);
        if (distance < min_distance) {
            min_distance = distance;
            closest_index = i;
        }
    }
    return closest_index;
}

uint32_t smallest_value_index(std::vector<double> values){
  double smallest_value = std::numeric_limits<double>::max();
  size_t index = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    if (values[i] < smallest_value) {
      smallest_value = values[i];
      index = i;
    }
  }
  return(index);
}

bool is_substring(const std::string& main_string, const std::string& sub_string) {
    // Find returns std::string::npos if the substring is not found
    return main_string.find(sub_string) != std::string::npos;
}

void deduplicate_solutions(std::vector<std::vector<uint32_t>> vectors, std::vector<std::vector<uint32_t>> &possible_permutations){
  for(uint32_t i=0; i < vectors.size(); i++){
    if(i == 0){
      possible_permutations.push_back(vectors[i]);
      continue;
    }
    bool add = true;
    for(uint32_t j=0; j < possible_permutations.size(); j++){
      bool same = std::equal(possible_permutations[j].begin(), possible_permutations[j].end(), vectors[i].begin());
      if(same && (possible_permutations[j].size() == vectors[i].size())){
        add = false;
        break;
      }
    }
    if(add){
      possible_permutations.push_back(vectors[i]);
    }
  }  
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
  double error = 0.05;
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

std::vector<std::vector<double>> solve_possible_solutions(std::vector<float> tmp_means, double error){
  uint32_t n = 6;
  std::vector<double> means;
  for(auto &x : tmp_means){
    for(uint32_t i=1; i < n; i++){
      if(x * i < 1+error && x > error){
        means.push_back(x);
      }
    }
    if(x <= error){
      means.push_back(x);
    }
  }
  std::vector<std::vector<double>> all_combinations;
  for (uint32_t i= 1; i <= means.size(); i++) {
    std::vector<double> current_combination;
    generate_combinations(means, current_combination, 0, i, all_combinations);
    
  }
  std::vector<std::vector<double>> viable_solutions;
  for(auto vec : all_combinations){
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0); 
    if(sum < 1+error && sum > 1-error && vec.size() > 1) {
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

void generate_permutations(const std::vector<uint32_t>& elements, int n, int target, std::vector<std::vector<uint32_t>> &other_tmp){
  std::vector<uint32_t> subset(elements);
  n = std::min(n, static_cast<int>(elements.size()));
  do {
    std::vector<uint32_t> tmp = std::vector<uint32_t>(subset.begin(), subset.begin() + n);
    int count = 0;
    for(auto val : tmp){
      if((int)val == target){
        count++;
      }
    }
    if(count >= 2){
      other_tmp.push_back(tmp);
    }
  } while (std::next_permutation(subset.begin(), subset.end()));
}

void noise_resampler(int n, uint32_t index, std::vector<std::vector<uint32_t>> &possible_permutations, uint32_t amount_resample){
  std::vector<uint32_t> tmp;
  for(uint32_t i=0; i < (uint32_t)n; i++){
    if(i == index){
      for(uint32_t j=0; j < amount_resample; j++){
        tmp.push_back(i);
      }
    } else {
      tmp.push_back(i);
    }
  }
  std::vector<std::vector<uint32_t>> other_tmp;
  generate_permutations(tmp, 2, index, other_tmp);
  generate_permutations(tmp, 3, index, other_tmp);
  generate_permutations(tmp, 4, index, other_tmp);
   
  deduplicate_solutions(other_tmp, possible_permutations);
  std::vector<std::vector<uint32_t>> keep_solutions;
  for(auto perm : possible_permutations){
    bool keep = false;
    for(auto clust : perm){
      if(clust != index) keep = true;
    }
    if(keep && perm.size() > 1){
      keep_solutions.push_back(perm);
    } else if(perm.size() == 1){
      keep_solutions.push_back(perm);
    }
  }
  possible_permutations.clear();
  possible_permutations = keep_solutions;
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
  double threshold = 2;
  std::vector<uint32_t> flagged_idx;
   
  for(uint32_t i=0; i < prob_matrix.size(); i++){
    double assigned_prob = prob_matrix[i][assigned[i]];
    std::vector<double> tmp = prob_matrix[i];
    tmp.erase(tmp.begin() + assigned[i]);
    std::sort(tmp.begin(), tmp.end(), std::greater<double>());
    for(uint32_t j=0; j < tmp.size(); j++){
      //std::cerr << "here " << tmp[j] << " " << assigned_prob << " " << exp(tmp[j]) << " " << exp(assigned_prob) << std::endl;
      //std::cerr << i << " " << j << std::endl;
      //if(assigned_prob < tmp[j]) continue; //indicates the assigned cluster was not the 'first' choice
      if(exp(tmp[j]) * threshold > exp(assigned_prob)){
        /*std::cerr << "HERE" << std::endl;
        std::cerr << tmp[j] << " assign " << assigned_prob << std::endl;
        std::cerr << exp(tmp[j]) << " exp assign" << exp(assigned_prob) << std::endl;*/
        flagged_idx.push_back(i);
      }
      break;
    }
    //std::cerr << "i bottom " << i << std::endl;
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
    /*for(auto p : permutations[i]){
      std::cerr << p << " ";
    }
    std::cerr << "\n";*/
    for(uint32_t j = 0; j < permutations[i].size(); j++){
      //num of variants in position must match permutation
      if(permutations[i].size() != prob_matrix.size()){ 
        score = -1000000;
        continue;
      }
      score += prob_matrix[j][permutations[i][j]];
    }
    //std::cerr << score << std::endl;
    scores.push_back(score);
  }
  uint32_t max_idx = std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));
  return(permutations[max_idx]);
}

void assign_variants_simple(std::vector<variant> &variants, std::vector<std::vector<double>> prob_matrix, uint32_t index, std::vector<double> means, uint32_t lower_n){
  double smallest_peak = *std::min_element(means.begin(), means.end()); 
  uint32_t n = prob_matrix.size();
  
  //find the unique positions
  std::vector<uint32_t> unique_pos;
  for(uint32_t i = 0; i < variants.size(); i++){
    if(variants[i].fake){
      std::vector<double> probs;
      for(uint32_t j=0; j < n; j++){
        probs.push_back(prob_matrix[j][i]);
      }
      auto max_it = std::max_element(probs.begin(), probs.end());
      uint32_t index = std::distance(probs.begin(), max_it);
      variants[i].cluster_assigned = index;
      continue;
    }
    if (std::find(unique_pos.begin(), unique_pos.end(), variants[i].position) == unique_pos.end()) {
      unique_pos.push_back(variants[i].position);
    }
  }
  //determine all possible permutations of assignments, disregard sum condition of E(u) ~= 1
  std::vector<std::vector<uint32_t>> possible_permutations;
  //what do we do here when we have fewer than two groups
  for(uint32_t i=1; i <= lower_n; i++){
    perm_generator(n, i, possible_permutations);
  }
  //if(smallest_peak < 0.15){
  noise_resampler(n, index, possible_permutations, 6);
  //}
 
  /*for(auto perm : possible_permutations){
    for(auto p : perm){
      std::cerr << p << " ";
    }
    std::cerr << "\n";
  }*/

  //now we loop every unique position and assign the max prob combo of variants
  for(uint32_t i=0; i < unique_pos.size(); i++){
    std::vector<uint32_t> pos_idxs;
    std::vector<std::vector<double>> tmp_prob;
    uint32_t j = 0;
    //all locations in the prob matrix for this position
    for(uint32_t k = 0;  k < variants.size(); k++){    
      if(variants[k].position  == unique_pos[i] && !variants[k].depth_flag){
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
    std::vector<uint32_t> assigned;
    assigned = calculate_joint_probabilities(tmp_prob, possible_permutations);
    if(assigned.size() == 0) continue;
    //make sure the assignment is concrete
    std::vector<uint32_t> assignment_flagged;
    assignment_flagged = compare_cluster_assignment(tmp_prob, assigned);

    for(uint32_t j=0; j < pos_idxs.size(); j++){
      std::vector<uint32_t>::iterator tmp = std::find(assignment_flagged.begin(), assignment_flagged.end(), j);
      uint32_t k = 0;
      for(uint32_t z =0; z < variants.size(); z++){
        //this pos was flagged as poorly assigned
        if(tmp != assignment_flagged.end() && k == pos_idxs[j]){
          //technically this could use work as it's repetitive
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
}

void determine_low_prob_positions(std::vector<variant> &variants){
  std::vector<uint32_t> low_prob_pos;
  for(uint32_t i=0; i< variants.size(); i++){
    if(variants[i].cluster_assigned == -1) continue;
    double prob = exp(variants[i].probabilities[variants[i].cluster_assigned]);
    double compare = (double)variants[i].depth / 1000000000;
    if(prob < compare){
      low_prob_pos.push_back(variants[i].position);
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
      } else if (std::abs(cm - variants[i].freq) > 0.10){
        //std::cerr << variants[i].freq << " " << " " << cm << std::endl;
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
      count++;
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
    count++;
    if(1/freq * depth < depth_cutoff || freq < lower_bound) continue;
    if (is_substring(nuc, "+") || is_substring(nuc, "-")) {
      deletion_positions.push_back(pos);
    }
    count++;
  }
  return(deletion_positions);
}

void parse_internal_variants(std::string filename, std::vector<variant> &variants, uint32_t depth_cutoff, float lower_bound, float upper_bound, std::vector<uint32_t> deletion_positions, uint32_t round_val){
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
    float freq = 0, gapped_freq=0, qual=0;
    std::string flag = "";
    std::string amp_flag = "";
    std::string primer_flag = "";
    std::string is_ref = "";
    split(line, '\t', row_values);
    pos = std::stoi(row_values[0]);
    bool del_pos = std::find(deletion_positions.begin(), deletion_positions.end(), pos) != deletion_positions.end();

    depth = std::stoi(row_values[2]);
    freq = std::stof(row_values[3]);
   
    float multiplier = pow(10, round_val);
    freq = round(freq * multiplier) / multiplier;
    
    gapped_freq = std::stof(row_values[4]);
    qual = std::stof(row_values[5]);
    flag = row_values[6];
    amp_flag = row_values[7];
    primer_flag = row_values[8];
    is_ref = row_values[9];
  
    variant tmp;
    tmp.position = pos;
    tmp.nuc = row_values[1];
    tmp.depth = depth;
    tmp.freq = freq;
    //tmp.transformed_freq = log(freq / (1 - freq));;
    tmp.transformed_freq = 1.0 / (1.0 + std::exp(-1 * (freq - 0.5)));
    tmp.qual = qual;
   
    if (flag == "TRUE"){
      tmp.amplicon_flux = true;
    } else {
      tmp.amplicon_flux = false;
    }
    if (amp_flag == "TRUE"){
      tmp.amplicon_masked = true;
    } else {
      tmp.amplicon_masked = false;
    }
    if (primer_flag == "TRUE"){
        tmp.primer_masked = false;
    } else {
        tmp.primer_masked = false;
    }
    if(1/gapped_freq * depth < depth_cutoff){
      if(pos == 3735){
        std::cerr << 1/freq * depth << std::endl;
      }
      tmp.depth_flag = true;
    } else {
      tmp.depth_flag = false;
    } 
    if (freq <= lower_bound || freq >= upper_bound){
      tmp.outside_freq_range = true;
    } else {
      tmp.outside_freq_range = false;
    }
    if (is_ref == "TRUE"){
      tmp.is_ref = true;
    } else {
      tmp.is_ref = false;
    }
    //TESTLINES
    if(qual < 20){
      tmp.qual_flag = true;
    } else {
      tmp.qual_flag = false;
    }
    if(del_pos && row_values[1] == "-" || row_values[1].find('+') != std::string::npos){
      tmp.del_flag = true;
    } else {
      tmp.del_flag = false;
    } 
    if(del_pos) {
      tmp.pos_del_flag = true;
    }
    variants.push_back(tmp); 
    count += 1;
  } 
}

std::vector<variant> gmm_model(std::string prefix, std::string output_prefix){
  uint32_t n=6;
  //changes based on error rate
  float lower_bound = 0.01;
  float upper_bound = 0.99;
  uint32_t depth_cutoff = 10;
  float quality_threshold = 20;
  uint32_t round_val = 4; 
  double universal_cluster = 0.97;
  double noise_cluster = 0.03;

  bool development_mode=true;

  float error_rate = cluster_error(prefix);
  lower_bound = error_rate;
  upper_bound = 1 - error_rate;
  std::cerr << lower_bound << " " << upper_bound << std::endl;

  //TESTLINES HEFTS CODE
  std::vector<std::string> heft_strings;
  std::vector<variant> base_variants;
  std::vector<uint32_t> deletion_positions = find_deletion_positions(prefix, depth_cutoff, lower_bound, upper_bound, round_val);

  parse_internal_variants(prefix, base_variants, depth_cutoff, lower_bound, upper_bound, deletion_positions, round_val);
  
  std::string filename = prefix + ".txt";

  //this whole things needs to be reconfigured
  uint32_t useful_var=0;
  std::vector<double> all_var;
  std::vector<variant> variants;
  std::vector<uint32_t> count_pos;

  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].del_flag && !base_variants[i].amplicon_masked && !base_variants[i].primer_masked){
      useful_var++;
      variants.push_back(base_variants[i]);
      all_var.push_back(base_variants[i].freq);      
      count_pos.push_back(base_variants[i].position);
      //std::cerr << base_variants[i].freq << " " << base_variants[i].position << " " << base_variants[i].nuc << std::endl;
    }
  }
  std::cerr << "useful var " << useful_var << std::endl;
  uint32_t lower_n = find_max_frequency_count(count_pos);
  std::cerr << "lower n " << lower_n << std::endl;
  if(useful_var == 0){
    variants.clear();
    return(variants);
  }
  //initialize armadillo dataset and populate with frequency data
  arma::mat data(1, useful_var, arma::fill::zeros);

  //(rows, cols) where each columns is a sample  
  uint32_t count=0;
  for(uint32_t i = 0; i < variants.size(); i++){   
    //check if variant should be filtered for first pass model
    double tmp = static_cast<double>(variants[i].freq); //transform
    data.col(count) = tmp;
    count += 1;
  }
  base_variants.clear();
  std::vector<double> means;
  std::vector<double> hefts;
  std::vector<uint32_t> exclude_cluster_indices;
  
  uint32_t counter = 2;
  uint32_t optimal_n = 0;
  gaussian_mixture_model retrained;

  //TESTLINES
  std::vector<std::vector<double>> mads;
  std::vector<std::vector<double>> percents;
  while(counter <= n){
    std::vector<double> tmp_mads;
    std::vector<double> tmp_percent_far;
    if(((float)useful_var / (float)counter) < 1){
      if(counter > 2){
        optimal_n = counter - 1;
      } else {
        optimal_n = counter;
      }
      break;
    }
    std::cerr << "\nn: " << counter << std::endl;
    retrained.means.clear();
    retrained.hefts.clear();
    retrained.prob_matrix.clear();
    retrained = retrain_model(counter, data, variants, lower_n, 0.001);
    bool optimal = true;
    assign_clusters(variants, retrained, lower_n);
    std::vector<std::vector<double>> clusters(counter);
    for(auto var : variants){
      int cluster = var.cluster_assigned;
      clusters[cluster].push_back(var.freq);
    }
    bool empty_cluster = false;
    for(auto data : clusters){
      double mean = calculate_mean(data);
      double mad = calculate_mad(data, mean);
      if(data.size() < 1){
        empty_cluster = true;
        std::cerr << "empty cluster" << std::endl;
        continue;
      }
      int count_far = 0;
      for(auto d : data){
        //std::cerr << d << " ";
        //this is for 2 standard dev
        if(std::abs(d-mean) > 0.05){
          count_far++;
        }
      }
      std::cerr << "\n";
      float percent_far = (float)count_far / (float)useful_var;
      std::cerr << "mean " << mean << " mad " << mad << " cluster size " << data.size() << " count far " << count_far << " percent far " << percent_far << std::endl;
       
      tmp_mads.push_back(mad);
      tmp_percent_far.push_back(percent_far);

      if(mad <= 0.10 && percent_far <= 0.10){
        optimal = true;
      } else {
        optimal = false;
        break;
      }
    }

    /*if(empty_cluster){
      optimal_n = counter-1;
      break;
    }*/

    if(optimal && !empty_cluster){
      optimal_n = counter;
      break;
    }
    counter++;     
    mads.push_back(tmp_mads);
    percents.push_back(tmp_percent_far);
  }
  std::cerr << "optimal n " << optimal_n << std::endl;
  retrained.means.clear();
  retrained.hefts.clear();
  retrained.prob_matrix.clear();
  retrained = retrain_model(optimal_n, data, variants, lower_n, 0.001);
  assign_clusters(variants, retrained, lower_n);
  std::vector<double> centroids;
  std::vector<std::vector<double>> c_clusters(optimal_n);
  for(auto var : variants){
    int cluster = var.cluster_assigned;
    c_clusters[cluster].push_back(var.freq);
  }
  for(auto data : c_clusters){
    double mean = calculate_mean(data);
    centroids.push_back(mean);
  }
  retrained = retrain_model_seeded(optimal_n, data, variants, lower_n, centroids);
  assign_clusters(variants, retrained, lower_n);
  std::vector<std::vector<double>> clusters(optimal_n);
  for(auto var : variants){
    int cluster = var.cluster_assigned;
    clusters[cluster].push_back(var.freq);
  }
  std::vector<double> final_means; 
  for(auto data : clusters){
    double mean = calculate_mean(data);
    /*for(auto d : data){
      std::cerr << d << " ";
    }
    std::cerr <<"\n";
    std::cerr << "end mean " << mean << std::endl;*/
    final_means.push_back(mean);
  }
  std::ofstream file;  
  //write mad to strings for use
  std::string mads_string;
  for(uint32_t l=0; l < mads.size(); l++){
    std::string tmp_str = "[";
    for(uint32_t d=0; d < mads[l].size(); d++){
      if(d != 0) tmp_str += ",";
      tmp_str += std::to_string(mads[l][d]);
    }
    tmp_str += "]\n";
    mads_string += tmp_str;
  }
  
  std::cerr << mads_string << std::endl;
  std::string mad_output = output_prefix + "_mad.txt";
  file.open(mad_output, std::ios::trunc);
  file << "MADS\n";
  file << mads_string;
  file.close(); 

  std::string percent_string;
  for(uint32_t l=0; l < percents.size(); l++){
    std::string tmp_str = "[";
    for(uint32_t d=0; d < percents[l].size(); d++){
      if(d != 0) tmp_str += ",";
      tmp_str += std::to_string(percents[l][d]);
    }
    tmp_str += "]\n";
    percent_string += tmp_str;
  }
  std::cerr << percent_string << std::endl;
  std::string percent_output = output_prefix + "_percent.txt";
  file.open(percent_output, std::ios::trunc);
  file << "PERCENTS\n";
  file << percent_string;
  file.close(); 

  //write means to string
  file.open(output_prefix + ".txt", std::ios::trunc);
  std::string means_string = "[";
  for(uint32_t j=0; j < final_means.size(); j++){
    if(j != 0) means_string += ",";
    means_string += std::to_string(final_means[j]);
  }
  means_string += "]";
  file << "means\n";
  file << means_string << "\n";
  file.close();

  //get the probability of each frequency being assigned to each gaussian
  //double prob_sum = 0;
  determine_low_prob_positions(variants);
  useful_var = 0;

  //look at all variants despite other parameters
  base_variants.clear();
  variants.clear();
 
  parse_internal_variants(prefix, base_variants, depth_cutoff, lower_bound, upper_bound, deletion_positions, round_val);
 
  calculate_gapped_frequency(base_variants, universal_cluster, noise_cluster);
 
  //could benefit from redoing lines in variants file as gapped/ungapped depth
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].outside_freq_range && !base_variants[i].pos_del_flag){
      useful_var += 1;
      variants.push_back(base_variants[i]);
    } else if(base_variants[i].pos_del_flag && !base_variants[i].gap_outside_freq_range){
      useful_var += 1;
      variants.push_back(base_variants[i]);
    }
  }
  //populate a new armadillo dataset with more frequencies
  arma::mat final_data(1, useful_var, arma::fill::zeros);
  count=0;
  for(uint32_t i = 0; i < variants.size(); i++){   
    double tmp;
    //check if variant should be filtered for first pass model
    if(!variants[i].pos_del_flag){
      tmp = static_cast<double>(variants[i].freq); //transform
    } else {
      tmp = static_cast<double>(variants[i].gapped_freq); //transform
    }
    final_data.col(count) = tmp;
    count += 1;
  }
  //recalculate the prob matrix based on the new dataset
  std::vector<std::vector<double>> prob_matrix;
  std::vector<double> tmp;
  for(uint32_t i=0; i < optimal_n; i++){
    arma::rowvec set_likelihood = retrained.model.log_p(final_data, i);
    tmp.clear();
    for(uint32_t j=0; j < final_data.n_cols; j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }
  retrained.prob_matrix = prob_matrix;
  assign_clusters(variants, retrained, lower_n);
  //determine_low_prob_positions(variants);
  //lets add back in the 100% variants
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].outside_freq_range && !base_variants[i].pos_del_flag){
      variants.push_back(base_variants[i]);
    } else if(base_variants[i].pos_del_flag && base_variants[i].gap_outside_freq_range){
      variants.push_back(base_variants[i]);
    }
  }
  if(development_mode){
    std::string cluster_output = output_prefix + "_cluster_data.txt";
    file.open(cluster_output, std::ios::trunc);
    file << "POS\tALLELE\tFREQ\tCLUSTER\tLIKELIHOOD\n";
    for(uint32_t i=0; i < variants.size(); i++){
      file << std::to_string(variants[i].position) << "\t"; 
      file << variants[i].nuc << "\t";
      if(!variants[i].pos_del_flag){
        file << std::to_string(variants[i].freq) << "\t";
      } else{
        file << std::to_string(variants[i].gapped_freq) << "\t";
      }
      file << std::to_string(variants[i].cluster_assigned) << "\t";
      if(variants[i].cluster_assigned != -1){
        float tmp = variants[i].probabilities[variants[i].cluster_assigned];
        file << std::to_string(tmp) << "\n"; 
      } else {
        file << "None\n";
      }
    }
    file.close();
  }
  return(variants);
}
