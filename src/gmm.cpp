#include "./include/armadillo"
#include "gmm.h"
#include "saga.h"
#include "call_consensus_clustering.h"
#include "solve_clustering.h"
#include "estimate_error.h"
#include "ref_seq.h"
#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>

std::vector<std::vector<double>> form_clusters(uint32_t n, std::vector<variant> variants){
  std::vector<std::vector<double>> clusters(n);
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].cluster_assigned != -1){
      clusters[variants[i].cluster_assigned].push_back(variants[i].gapped_freq);
    }   
  }
  return(clusters);
}

double calculate_mean(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0f;
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0f);
    return sum / data.size();
}

uint32_t find_max_frequency_count(const std::vector<uint32_t>& nums) {
    std::unordered_map<uint32_t, uint32_t> frequency_map;
    uint32_t max_count = 0;

    for (const auto& num : nums) {
        frequency_map[num]++;
        if (frequency_map[num] > max_count) {
            max_count = frequency_map[num];
        }
    }

    return max_count;
}

double calculate_mad(const std::vector<double>& data, double mean){
    double absDevSum = 0.0;
    for (double value : data) {
        absDevSum += std::abs(value - mean);
    }
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
  //handle everything but insertions
  assign_variants_simple(variants, gmodel.prob_matrix, index, lower_n, false);
  //handle insertions
  assign_variants_simple(variants, gmodel.prob_matrix, index, lower_n, true);
}

/**
* @brief Train a KMeans model to seed other analyses.
* @params n An integer indicating the number of model components.
* @param error A booleans value that indicates if this kmeans is being used to detected error levels.
* @return kmeans_model A kmeans_modle object storing centroids and clusters.
*/
kmeans_model train_model(uint32_t n, arma::mat data, bool error) {
  arma::mat centroids;
  arma::mat initial_means(1, n, arma::fill::zeros);
  kmeans_model model;

  std::vector<std::vector<double>> all_centroids;
  std::vector<double> total_distances; 
  std::vector<double> std_devs;
  for(uint32_t j=0; j < 15; j++){
    bool status2 = arma::kmeans(centroids, data, n, arma::random_spread, 10, false);
    if(!status2) continue;
    double total_dist = 0;
    for(auto point : data){
      //using std::min_element to find the closest element
      auto closest_it = std::min_element(centroids.begin(), centroids.end(),
        [point](double a, double b) {
            return std::abs(a - point) < std::abs(b - point);
        });
       uint32_t index = std::distance(centroids.begin(), closest_it);      
       total_dist += std::abs(point-centroids[index]);
    }  
    std::vector<double> tmp = arma::conv_to<std::vector<double>>::from(centroids);
    double stddev = calculate_standard_deviation(tmp);
    std_devs.emplace_back(stddev);
    all_centroids.emplace_back(std::move(tmp));
    total_distances.emplace_back(total_dist);
  }
  uint32_t index;
  if(!error){
    auto min_it = std::min_element(total_distances.begin(), total_distances.end());
    index = std::distance(total_distances.begin(), min_it);
  } else {
    auto max_it = std::max_element(std_devs.begin(), std_devs.end());
    index = std::distance(std_devs.begin(), max_it);
  }
  std::vector<double> centroid_vec = all_centroids[index];

  //assign points to cluster
  std::vector<std::vector<double>> clusters(n);
  for(auto point : data){
    //using std::min_element to find the closest element
    auto closest_it = std::min_element(centroid_vec.begin(), centroid_vec.end(),
        [point](double a, double b) {
            return std::abs(a - point) < std::abs(b - point);
        });
    uint32_t index = std::distance(centroid_vec.begin(), closest_it);
    clusters[index].push_back(point);  
  }
  model.n = n;
  model.means = centroid_vec;
  model.clusters = clusters;
  return(model);
}

//function used for production
gaussian_mixture_model retrain_model(uint32_t n, arma::mat data, std::vector<variant> &variants, uint32_t lower_n, double var_floor){
  double initial_covariance = 0.005;
  gaussian_mixture_model gmodel; 
  gmodel.n = n;
  arma::mat initial_means(1, n, arma::fill::zeros);

  arma::mat cov (1, n, arma::fill::zeros);
  std::vector<double> total_distances;
  std::vector<std::vector<double>> all_centroids;
  
  //run a kmeans to seed the GMM
  kmeans_model initial_model = train_model(n, data, false);

  for(uint32_t c=0; c < initial_model.means.size(); c++){
    initial_means.col(c) = (double)initial_model.means[c];
    cov.col(c) = initial_covariance;
  }
  
  arma::gmm_diag model;
  model.reset(1, n);
  model.set_means(initial_means);
  model.set_dcovs(cov);
  bool status = model.learn(data, n, arma::eucl_dist, arma::keep_existing, 1, 10, var_floor, false);
  if(!status){
    std::cerr << "GMM failed to converge" << std::endl;
    exit(1);
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
  assign_clusters(variants, gmodel, lower_n);

  std::vector<std::vector<double>> clusters = form_clusters(n, variants);
  gmodel.clusters = clusters;
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
        if(std::abs(theoretical_peak - mean) < error){
          useful = true;
          break;
        }  
      }
      if(!useful){
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

void generate_permutations(const std::vector<uint32_t>& elements, int n, int target, std::vector<std::vector<uint32_t>> &other_tmp) {
    std::vector<uint32_t> subset(elements);
    n = std::min(n, static_cast<int>(elements.size()));
    std::set<std::vector<uint32_t>> seen; // store unique permutations
    std::sort(subset.begin(), subset.end()); // needed for std::next_permutation
    do {
        std::vector<uint32_t> tmp(subset.begin(), subset.begin() + n);
        int count = std::count(tmp.begin(), tmp.end(), target);
        if (count >= 2 && seen.insert(tmp).second) {
            // only add if this permutation is new
            other_tmp.push_back(tmp);
        }
    } while (std::next_permutation(subset.begin(), subset.end()));
}

void noise_resampler(uint32_t n, uint32_t index, std::vector<std::vector<uint32_t>> &possible_permutations, uint32_t amount_resample){
  std::vector<uint32_t> tmp;
  for(uint32_t i=0; i < n; i++){
    if(i == index){
      for(uint32_t j=0; j < amount_resample; j++){
        tmp.push_back(i);
      }
    } else {
      tmp.push_back(i);
    }
  }
  generate_permutations(tmp, 2, index, possible_permutations);
  generate_permutations(tmp, 3, index, possible_permutations);
  generate_permutations(tmp, 4, index, possible_permutations);
  
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
      if(exp(tmp[j]) * threshold > exp(assigned_prob)){
        flagged_idx.push_back(i);
      }
      break;
    }
  }
  return(flagged_idx);
}



/**
 * @brief Selects the permutation of assignments that maximizes the joint probability.
 *
 * This function evaluates a set of possible assignments (permutations) and computes the
 * total joint probability score for each, using the provided probability matrix.
 * It returns the permutation that yields the highest score.
 *
 * @param prob_matrix A 2D vector of probabilities, sized [n_variants][n_clusters].
 * @param permutations A vector of permutations to evaluate, each representing a possible assignment.
 *                     Each permutation must have a size equal to the number of clusters.
 * @return The permutation (as a vector of cluster indices) with the highest joint probability.
 */
std::vector<uint32_t> calculate_joint_probabilities(const std::vector<std::vector<double>> prob_matrix, const std::vector<std::vector<uint32_t>> permutations) {
  if (permutations.empty() || prob_matrix.empty()) {
    return {};
  }
  size_t n_clusters = prob_matrix.size();
  double best_score = -std::numeric_limits<double>::infinity(); 
  size_t best_index = 0;
  for (size_t i = 0; i < permutations.size(); ++i) {
    const auto& perm = permutations[i];
    if (perm.size() != n_clusters) {
        continue;
    }
    double score = 0.0;
    for (size_t j = 0; j < n_clusters; ++j) {
      // Guard against invalid index in permutation
      if (perm[j] >= prob_matrix[j].size()) {
        score = -std::numeric_limits<double>::infinity();
        break;
      }
      score += prob_matrix[j][perm[j]];
    }
    if (score > best_score) {
      best_score = score;
      best_index = i;
    }
  }
  return permutations[best_index];
}

void assign_variants_simple(std::vector<variant> &variants, std::vector<std::vector<double>> prob_matrix, uint32_t index, uint32_t lower_n, bool insertions) {
  uint32_t n = prob_matrix.size();
  std::unordered_map<uint32_t, std::vector<std::string>> all_nts;
  std::unordered_map<uint32_t, std::vector<uint32_t>> pos_to_variant_indices;

  // Map positions to variant indices and nucleotides
  for (uint32_t i = 0; i < variants.size(); ++i) {
    uint32_t pos = variants[i].position;
    all_nts[pos].push_back(variants[i].nuc);
    pos_to_variant_indices[pos].push_back(i);
  }

  std::vector<uint32_t> unique_pos;
  for (const auto& kv : all_nts)
    unique_pos.push_back(kv.first);

  // Generate all permutations up to lower_n
  std::vector<std::vector<uint32_t>> possible_permutations;
  for (uint32_t i = 1; i <= lower_n; ++i)
    perm_generator(n, i, possible_permutations);

  noise_resampler(n, index, possible_permutations, 6);

  // Assignment by position
  for (uint32_t pos : unique_pos) {
    std::vector<uint32_t> pos_idxs;
    std::vector<std::vector<double>> tmp_prob;

    for (uint32_t variant_idx : pos_to_variant_indices[pos]) {
      auto& var = variants[variant_idx];
      if ((var.nuc.find('+') != std::string::npos && !insertions)|| var.depth_flag)
        continue;
      else if ((var.nuc.find('+') == std::string::npos && insertions)|| var.depth_flag)
        continue;
      pos_idxs.push_back(variant_idx);
      std::vector<double> prob_column;
      prob_column.reserve(n);
      for (uint32_t row = 0; row < n; ++row)
        prob_column.push_back(prob_matrix[row][variant_idx]);

      tmp_prob.push_back(std::move(prob_column));
    }

    if (pos_idxs.empty())
      continue;

    std::vector<uint32_t> assigned = calculate_joint_probabilities(tmp_prob, possible_permutations);
    if (assigned.empty())
      continue;

    std::vector<uint32_t> assignment_flagged = compare_cluster_assignment(tmp_prob, assigned);

    for (uint32_t i = 0; i < pos_idxs.size(); ++i) {
      uint32_t v_idx = pos_idxs[i];
      if (std::find(assignment_flagged.begin(), assignment_flagged.end(), i) != assignment_flagged.end()) {
        variants[v_idx].vague_assignment = true;
      }
      variants[v_idx].cluster_assigned = assigned[i];
    }
  }
}

void split(std::string &s, char delim, std::vector<std::string> &elems){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<double> split_csv_double(const std::string& input) {
    std::vector<double> result;
    std::stringstream ss(input);
    std::string token;
    
    while (std::getline(ss, token, ',')) {
      result.push_back(std::stod(token));
    }
    return result;
}


std::vector<uint32_t> split_csv(const std::string& input) {
    std::vector<uint32_t> result;
    std::stringstream ss(input);
    std::string token;
    
    while (std::getline(ss, token, ',')) {
      result.push_back(std::stoi(token));
    }
    return result;
}

void set_freq_range_flags(std::vector<variant> &variants, double lower_bound, double upper_bound){
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].gapped_freq <= lower_bound || variants[i].gapped_freq >= upper_bound){
      variants[i].outside_freq_range = true;
    } else {
      variants[i].outside_freq_range = false;
    }
  }
}

void deletion_variant(std::vector<std::string> row_values, std::vector<variant> &variants, uint32_t depth_cutoff, double compare_quality, double multiplier){
  std::string deletion = row_values[3];
  for(uint32_t i=1; i < deletion.size(); i++){
    variant tmp;
    tmp.position = std::stoi(row_values[1]) + i; 
    tmp.nuc = '-';
    tmp.depth = std::stoi(row_values[7]);
    tmp.total_depth = std::stoi(row_values[11]);
    double freq = std::stod(row_values[10]);
    tmp.freq = std::round(freq * multiplier) / multiplier;
    tmp.qual = std::stod(row_values[6]);
    auto to_bool = [](const std::string& s) -> bool {return s == "TRUE" || s == "true" || s == "1";};
    if(row_values.size() > 20){
      tmp.gapped_depth = std::stoi(row_values[21]);
      tmp.gapped_freq = std::round(std::stod(row_values[20]) * multiplier) / multiplier;
      tmp.amplicon_flux = to_bool(row_values[22]);
      tmp.amplicon_masked = to_bool(row_values[23]);
      tmp.primer_masked = to_bool(row_values[24]);
      tmp.std_dev = std::stod(row_values[25]);
      tmp.amplicon_numbers = split_csv(row_values[27]);
      tmp.freq_numbers = split_csv_double(row_values[26]);
      tmp.version_1_var = false;
    } else {
      tmp.gapped_freq = 0.0;
      tmp.amplicon_flux = false;
      tmp.amplicon_masked = false;
      tmp.primer_masked = false;
      tmp.std_dev = 0;
      tmp.version_1_var = true;
    }
    if(!tmp.version_1_var && tmp.total_depth < depth_cutoff){
      tmp.depth_flag = true;
    } else if(tmp.version_1_var && tmp.depth < depth_cutoff){
      tmp.depth_flag = true;
    } else {
      tmp.depth_flag = false;
    }
    if(tmp.qual < compare_quality){
      tmp.qual_flag = true;
    }   else {
      tmp.qual_flag = false;
    }
    variants.push_back(std::move(tmp)); 
  }
}

/**
 * @brief Parses an internally formatted variant file and populates a vector of variant objects.
 *
 * This function reads a tab-delimited file line-by-line (skipping the header),
 * extracts relevant variant information, and fills a vector of `variant` structs.
 * It supports both version 1 and newer variant file formats with additional metadata.
 *
 * @param filename           Path to the tab-delimited variant file.
 * @param variants           Reference to a vector where parsed variant entries will be stored.
 * @param depth_cutoff       Minimum depth threshold (currently unused in this function).
 * @param round_val          Number of decimal places to round frequencies to.
 * @param quality_threshold  Minimum quality score threshold (currently unused in this function).
 *
 * @note The function assumes a specific column structure in the input file.
 *       It handles both older (≤20 columns) and newer (>20 columns) variant formats.
 *       Insertions and deletions may be represented with special notations in the input.
 *
 * @warning Malformed or incomplete lines (with <12 columns) are silently skipped.
 *
 * @see variant
 */
void parse_internal_variants(std::string filename, std::vector<variant> &variants, uint32_t depth_cutoff, uint32_t round_val, uint8_t quality_threshold, std::string reference_file){
  std::ifstream infile(filename);
  std::string line;
  uint32_t count = 0;
  double multiplier = pow(10, round_val);
  double compare_quality = static_cast<double>(quality_threshold);

  auto to_bool = [](const std::string& s) -> bool {return s == "TRUE" || s == "true" || s == "1";};
  //track which ref alleles we've already added
  //std::vector<uint32_t> ref_pos_used;
  while (std::getline(infile, line)) {
    if(count++ == 0) continue;
    std::vector<std::string> row_values;
    split(line, '\t', row_values);

    variant tmp;
    tmp.position = std::stoi(row_values[1]);
    tmp.nuc = row_values[3];
    bool is_ins = std::find(tmp.nuc.begin(), tmp.nuc.end(), '+') != tmp.nuc.end();
    bool is_del = std::find(tmp.nuc.begin(), tmp.nuc.end(), '-') != tmp.nuc.end();
    if(is_ins || is_del) {
      tmp.position = tmp.position + 1;
    }
    //seperate deletions into individual reference based variants
    if(is_del){
      deletion_variant(row_values, variants, depth_cutoff, compare_quality, multiplier);
      continue;
    }
    tmp.depth = std::stoi(row_values[7]);
    tmp.total_depth = std::stoi(row_values[11]);
    tmp.freq = std::round(std::stof(row_values[10]) * multiplier) / multiplier;
    tmp.qual = std::stod(row_values[9]);

    if(row_values.size() > 20){
      tmp.gapped_freq = round(std::stod(row_values[20]) * multiplier) / multiplier;
      tmp.gapped_depth = std::stoi(row_values[21]);
      tmp.amplicon_flux = to_bool(row_values[22]);
      tmp.amplicon_masked = to_bool(row_values[23]);
      tmp.primer_masked = to_bool(row_values[24]);
      tmp.std_dev = std::stod(row_values[25]);
      tmp.amplicon_numbers = split_csv(row_values[27]);
      tmp.freq_numbers = split_csv_double(row_values[26]);
      tmp.version_1_var = false;
    } else {
      tmp.gapped_freq = 0.0;
      tmp.amplicon_flux = false;
      tmp.amplicon_masked = false;
      tmp.primer_masked = false;
      tmp.std_dev = 0;
      tmp.version_1_var = true;
    }

    if(tmp.total_depth < depth_cutoff){
      tmp.depth_flag = true;
    } else {
      tmp.depth_flag = false;
    }
    if(tmp.qual < compare_quality){
      tmp.qual_flag = true;
    } else {
      tmp.qual_flag = false;
    }
    variants.push_back(std::move(tmp)); 
  } 
}

std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, std::vector<double> &solution, std::vector<double> &means, std::string ref){
  if(ref.empty()){
    std::cerr << "Please provide a reference sequence." << std::endl;
    exit(1);
  }
  uint32_t n=8;
  uint32_t round_val = 4; 
  bool development_mode=true;

  std::vector<variant> base_variants;
  parse_internal_variants(prefix, base_variants, min_depth, round_val, min_qual, ref);
  double error_rate = cluster_error(base_variants, min_qual, min_depth);
  double lower_bound = 1-error_rate+0.0001;
  double upper_bound = error_rate-0.0001;
  set_freq_range_flags(base_variants, lower_bound, upper_bound);
  std::cerr << "gmm model lower " << lower_bound << " upper " << upper_bound << std::endl;
   
  //this whole things needs to be reconfigured
  uint32_t useful_var=0;
  std::vector<variant> variants;
  std::vector<uint32_t> count_pos;

  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked){
      useful_var++;
      variants.push_back(base_variants[i]);
      count_pos.push_back(base_variants[i].position);
      //std::cerr << base_variants[i].freq << " " << base_variants[i].position << " " << base_variants[i].nuc << " " << base_variants[i].depth << " " << base_variants[i].gapped_freq << std::endl;
    }
  }
  std::cerr << "useful var " << useful_var << std::endl;

  if(useful_var < 2){
    variants.clear();
    return(variants);
  }
  uint32_t lower_n = find_max_frequency_count(count_pos);
  //initialize armadillo dataset and populate with frequency data
  arma::mat data(1, useful_var, arma::fill::zeros);

  //(rows, cols) where each columns is a sample  
  uint32_t count=0;
  for(uint32_t i = 0; i < variants.size(); i++){   
    double tmp = static_cast<double>(variants[i].gapped_freq);
    data.col(count) = tmp;
    count += 1;
  }
  std::vector<uint32_t> exclude_cluster_indices;
  
  uint32_t counter = 2;
  uint32_t optimal_n = 0;
  gaussian_mixture_model retrained;

  //TESTLINES
  while(counter <= n){
    if(((double)useful_var / (double)counter) < 1){
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
    std::vector<std::vector<double>> clusters = form_clusters(counter, variants);
    bool empty_cluster = false;
    for(auto data : clusters){
      double mean = calculate_mean(data);
      double mad = calculate_mad(data, mean);
      if(data.size() < 1){
        empty_cluster = true;
        std::cerr << "empty cluster" << std::endl;
        continue;
      }
      std::cerr << "mean " << mean << " mad " << mad << " cluster size " << data.size() << " ratio " << std::endl;
      if(data.size() > 5){
        if(mad <= 0.10){
          optimal = true;
        } else {
          optimal = false;
          break;
        }
      } else {
        if(mad <= 0.03){
          optimal = true;
        }  else{
          optimal = false;
          break;
        }
      }
    }

    if(optimal && !empty_cluster){
      optimal_n = counter;
      break;
    }
    counter++;     
  }
  if(optimal_n != retrained.means.size()){
    retrained.means.clear();
    retrained.hefts.clear();
    retrained.prob_matrix.clear();
    retrained = retrain_model(optimal_n, data, variants, lower_n, 0.001);
  }
  for(auto cluster : retrained.clusters){
    double mean = calculate_mean(cluster);
    means.push_back(mean);
  }
  retrained.means = means;
  //TODO this can be put in function
  std::ofstream file;
  if(development_mode){
    std::vector<std::vector<double>> clusters = form_clusters(optimal_n, variants);
    std::vector<double> final_means; 
    for(auto data : clusters){
      double mean = calculate_mean(data);
      final_means.push_back(mean);
    }
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
  }

  //get the probability of each frequency being assigned to each gaussian
  useful_var = 0;

  //look at all variants despite other parameters
  variants.clear();
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].outside_freq_range){
      useful_var += 1;
      variants.push_back(base_variants[i]);
    }
  }
  //populate a new armadillo dataset with more frequencies
  arma::mat final_data(1, useful_var, arma::fill::zeros);
  for(uint32_t i = 0; i < variants.size(); i++){   
    final_data.col(i) = static_cast<double>(variants[i].gapped_freq);
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
 
  //lets add back in the 100% variants
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].outside_freq_range){
      variants.push_back(base_variants[i]);
    }
  }
  solve_clusters(variants, retrained, error_rate, solution);
  return(variants);
}
