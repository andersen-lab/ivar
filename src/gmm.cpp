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

uint32_t elbow_method(std::vector<double> ics,
                      std::vector<double> ns,
                      std::vector<double> exclude_ns) {
    std::vector<double> slopes;
    std::vector<uint32_t> valid_ns;

    // Build aligned slope and n vectors, skipping excluded ns
    for (size_t i = 0; i + 1 < ns.size(); ++i) {
        double slope = (ics[i+1] - ics[i]) / (ns[i+1] - ns[i]);
        uint32_t n_val = static_cast<uint32_t>(ns[i+1]);
          slopes.push_back(slope);
          valid_ns.push_back(n_val);
    }

    if (slopes.empty()) {
        std::cerr << "Warning: no valid slopes after exclusion, defaulting to n=1" << std::endl;
        return 1;
    }

    // Find largest increase (smallest -> largest slope transition)
    double max_jump = -1e9;
    uint32_t jump_index = 0;
    for (uint32_t i = 0; i + 1 < slopes.size(); ++i) {
        double jump = slopes[i+1] - slopes[i];
        if (jump > max_jump) {
            max_jump = jump;
            jump_index = i;  // return the "smaller" slope's index
        }
    }

    std::cerr << "Elbow at n=" << valid_ns[jump_index]
              << " with slope change " << slopes[jump_index]
              << " -> " << slopes[jump_index+1] << std::endl;

    return valid_ns[jump_index];
}

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

void generate_ordered(const std::vector<uint32_t>& elements,
                      uint32_t n,
                      uint32_t target,
                      std::vector<std::vector<uint32_t>> &results) {
    std::vector<uint32_t> seq;
    seq.reserve(n);
    int targetCount = 0;

    std::function<void()> backtrack = [&]() {
        if (seq.size() == n) {
            if (targetCount >= 2) results.push_back(seq);
            return;
        }
        for (uint32_t e : elements) {
            seq.push_back(e);
            if (e == target) targetCount++;
            backtrack();
            if (e == target) targetCount--;
            seq.pop_back();
        }
    };

    backtrack();
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
std::vector<uint32_t> calculate_joint_probabilities(const std::vector<std::vector<double>>  &prob_matrix, const std::vector<std::vector<uint32_t>> &permutations) {
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

void assign_variants_simple(std::vector<variant> &variants, 
                            const std::vector<std::vector<double>> &prob_matrix, 
                            bool insertions, 
                            bool &clustering_failed, 
                            const std::vector<std::vector<uint32_t>> &possible_permutations,
                            const std::vector<uint32_t> &unique_pos, 
                            std::unordered_map<uint32_t, std::vector<uint32_t>> &pos_to_variant_indices) {
  uint32_t n = prob_matrix.size();


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

    //here we have more variants trying to assign than clusters
    if(tmp_prob.size() > assigned.size()) {
      std::cerr << "HERE " << tmp_prob.size() << " " << assigned.size() << " " << pos << std::endl;
      clustering_failed = true;
      return;
    }

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

void assign_clusters(std::vector<variant> &variants, 
                    gaussian_mixture_model gmodel, 
                    bool &clustering_failed, 
                    std::vector<std::vector<uint32_t>> possible_permutations,
                    std::vector<uint32_t> unique_pos,
                    std::unordered_map<uint32_t, std::vector<uint32_t>> pos_to_variant_indices
                    ){
  std::vector<std::vector<double>> tv = transpose_vector(gmodel.prob_matrix);
  uint32_t j = 0;
  for(uint32_t i=0; i < variants.size(); i++){
    variants[i].probabilities = tv[j];
    j++;
  }
  uint32_t index = smallest_value_index(gmodel.means);
  //generate all permutations up to lower_n
  noise_resampler(gmodel.n, index, possible_permutations, 6);
  
  //handle everything but insertions
  assign_variants_simple(variants, gmodel.prob_matrix, false, clustering_failed, possible_permutations, unique_pos, pos_to_variant_indices);
  //handle insertions
  assign_variants_simple(variants, gmodel.prob_matrix, true, clustering_failed, possible_permutations, unique_pos, pos_to_variant_indices);
}

void assign_all_variants(std::vector<variant> &variants, 
                        std::vector<variant> base_variants, 
                        gaussian_mixture_model &gmodel, 
                        double lower_bound, 
                        double upper_bound) {
  std::vector<variant> tmp_var;
  uint32_t count = 0;

  for(uint32_t i = 0; i < base_variants.size(); i++){
    //previously not assigned due to possible amplicon flux
    if(base_variants[i].qual_flag) continue;
    if((!base_variants[i].outside_freq_range && (base_variants[i].amplicon_flux || base_variants[i].amplicon_masked || !base_variants[i].include_clustering)) || (base_variants[i].outside_freq_range && base_variants[i].gapped_freq >= lower_bound && base_variants[i].gapped_freq <= upper_bound)){
      count++;
      tmp_var.push_back(base_variants[i]);
      variants.push_back(base_variants[i]);
    }
  }
  std::unordered_map<uint32_t, std::vector<std::string>> all_nts;
  std::unordered_map<uint32_t, std::vector<uint32_t>> pos_to_variant_indices;
  for (uint32_t i = 0; i < variants.size(); ++i) {
    uint32_t pos = variants[i].position;
    all_nts[pos].push_back(variants[i].nuc);
    pos_to_variant_indices[pos].push_back(i);
  }

  std::vector<uint32_t> unique_pos;
  for (const auto& kv : all_nts)
    unique_pos.push_back(kv.first);

  //populate a new armadillo dataset with more frequencies
  arma::mat final_data(1, count, arma::fill::zeros);
  for(uint32_t i = 0; i < tmp_var.size(); i++){
    final_data.col(i) = static_cast<double>(tmp_var[i].gapped_freq);
  }
  //recalculate the prob matrix based on the new dataset
  std::vector<std::vector<double>> prob_matrix;
  std::vector<double> tmp;
  for(uint32_t i=0; i < gmodel.n; i++){
    arma::rowvec set_likelihood = gmodel.model.log_p(final_data, i);
    tmp.clear();
    for(uint32_t j=0; j < final_data.n_cols; j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }
  std::vector<std::vector<double>> stacked_matrix;
  for (size_t i = 0; i < gmodel.prob_matrix.size(); ++i) {
    std::vector<double> row = gmodel.prob_matrix[i];
    row.insert(row.end(), prob_matrix[i].begin(), prob_matrix[i].end());
    stacked_matrix.push_back(row);
  }
  gmodel.prob_matrix = stacked_matrix;
  bool clustering_failed = false;
  std::vector<std::vector<uint32_t>> possible_permutations;
  for (uint32_t i = 1; i <= gmodel.lower_n; ++i) {
    perm_generator(gmodel.n, i, possible_permutations);
  }
  assign_clusters(variants, gmodel, clustering_failed, possible_permutations, unique_pos, pos_to_variant_indices);
  if(clustering_failed){
    //std::cerr << "clustering failed" << std::endl;
    //exit(1);
  }
}

void add_noise_variants(std::vector<variant> &variants, std::vector<variant> base_variants){
  //lets add back in the 100% variants
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].outside_freq_range){
      variants.push_back(base_variants[i]);
    }
  }
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
 
  std::vector<double> means;
  std::vector<std::vector<double>> clusters(n);
  bool status = true;
  status = arma::kmeans(centroids, data, n, arma::random_subset, 10, false);
  if(!status) return(model);
  for(uint32_t c=0; c < centroids.n_cols; c++){
    means.push_back(centroids(0,c));
  }
  model.n = n;
  model.means = means;
  return(model);
}

gaussian_mixture_model retrain_model(uint32_t n, 
                                    arma::mat data, 
                                    arma::mat second_dim_data, 
                                    std::vector<variant> &variants, 
                                    uint32_t lower_n, 
                                    double var_floor, 
                                    bool &clustering_failed){

   //this is used in the variant assignement portion of the code
  std::unordered_map<uint32_t, std::vector<std::string>> all_nts;
  std::unordered_map<uint32_t, std::vector<uint32_t>> pos_to_variant_indices;
  //map positions to variant indices and nucleotides
  for (uint32_t i = 0; i < variants.size(); ++i) {
    uint32_t pos = variants[i].position;
    all_nts[pos].push_back(variants[i].nuc);
    pos_to_variant_indices[pos].push_back(i);
  }

  std::vector<uint32_t> unique_pos;
  for (const auto& kv : all_nts)
    unique_pos.push_back(kv.first);
  
  double initial_covariance = 0.005;
  uint32_t j = 0;
  std::vector<gaussian_mixture_model> all_models;
  std::vector<double> all_bics;

  while(j < 10) {
    gaussian_mixture_model gmodel;
    gmodel.n = n;
    gmodel.lower_n = lower_n;
    arma::mat initial_means(1, n, arma::fill::zeros);

    arma::mat cov (1, n, arma::fill::zeros);
    std::vector<double> total_distances;
    std::vector<std::vector<double>> all_centroids;

    //run a kmeans to seed the GMM
    kmeans_model initial_model = train_model(n, second_dim_data, false);

    for(uint32_t c=0; c < n; c++){
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
      clustering_failed = true;
      continue;
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
      for(uint32_t k=0; k < data.n_cols; k++){
        tmp.push_back((double)set_likelihood[k]);
      }
      prob_matrix.push_back(tmp);
    }
    gmodel.dcovs = dcovs;
    gmodel.prob_matrix = prob_matrix;
    gmodel.means = means;
    gmodel.hefts = hefts;
    gmodel.model = model;
    std::vector<std::vector<uint32_t>> possible_permutations;
    for (uint32_t i = 1; i <= gmodel.lower_n; ++i){
      perm_generator(gmodel.n, i, possible_permutations);
    }
    assign_clusters(variants, gmodel, clustering_failed, possible_permutations, unique_pos, pos_to_variant_indices);

    std::vector<std::vector<double>> clusters = form_clusters(n, variants);
    gmodel.clusters = clusters;

    //cluster probabilities
    std::vector<std::vector<double>> cluster_prob(n);
    for(auto var : variants){
      uint32_t assigned = var.cluster_assigned;
      if(assigned != (uint32_t)-1){
        double prob = var.probabilities[assigned];
        cluster_prob[assigned].push_back(prob);
      }
    }
    double total_likelihoods = 0;
    for(auto val : cluster_prob){
      double sum = std::accumulate(val.begin(), val.end(), 0.0);
      total_likelihoods += sum;
      gmodel.cluster_probabilities.push_back(sum/(double)val.size());
    }
    double bic = ((2 * n) + (n-1)) * std::log(data.size()) - (2 * total_likelihoods);
    gmodel.bic = bic;

    all_bics.push_back(bic);
    all_models.push_back(gmodel);
    j++;
  }
  auto it = std::min_element(all_bics.begin(), all_bics.end());
  uint32_t index_best = std::distance(all_bics.begin(), it);
  gaussian_mixture_model chosen = all_models[index_best];

  return(chosen);
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

void noise_resampler(uint32_t n, uint32_t index, std::vector<std::vector<uint32_t>> &possible_permutations, uint32_t amount_resample) {
  std::vector<uint32_t> tmp;
  tmp.reserve(n + amount_resample);

  for (uint32_t i = 0; i < n; i++) {
    if (i == index)
      tmp.insert(tmp.end(), amount_resample, i);
      else
        tmp.push_back(i);
  }

  generate_ordered(tmp, 2, index, possible_permutations);
  generate_ordered(tmp, 3, index, possible_permutations);
  generate_ordered(tmp, 4, index, possible_permutations);


  possible_permutations.erase(
    std::remove_if(possible_permutations.begin(),
                    possible_permutations.end(), [&](const std::vector<uint32_t>& perm) {
                    if (perm.size() == 1) return false;
                    return std::all_of(perm.begin(), perm.end(),
                    [&](uint32_t v){ return v == index; });
    }),
  possible_permutations.end());
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

void set_advanced_freq_range_flags(std::vector<variant> &variants, uint32_t max_pos, double lower_bound) {
  std::unordered_map<uint32_t, std::vector<variant>> variants_by_pos;
  for (auto &var : variants) {
    variants_by_pos[var.position].push_back(var);
  }

  for (auto &[pos, vec] : variants_by_pos) {
    if (pos >= max_pos) continue;

    double total_noise = 0.0;
    for (const auto &v : vec) {
      if (v.outside_freq_range) {
        total_noise += v.gapped_freq;
      }
    }

    if (total_noise < lower_bound) continue;

    for (auto &v : vec) {
      if (!v.outside_freq_range && (1.0 - v.gapped_freq) >= total_noise && v.gapped_freq > 0.5) {
        v.outside_freq_range = true;
      }
    }
  }
}


void set_freq_range_flags(std::vector<variant> &variants, double lower_bound, double upper_bound, bool advanced){
  uint32_t max_pos = 0;
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].position > max_pos) max_pos = variants[i].position;
    if(variants[i].gapped_freq <= lower_bound || variants[i].gapped_freq >= upper_bound){
      variants[i].outside_freq_range = true;
    } else {
      variants[i].outside_freq_range = false;
    }
  }
  if(advanced){
    set_advanced_freq_range_flags(variants, max_pos, lower_bound);
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

void parse_internal_variants(std::string filename, std::vector<variant> &variants, uint32_t depth_cutoff, uint32_t round_val, uint8_t quality_threshold){
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
    tmp.nuc = row_values[3];
    tmp.position = std::stoi(row_values[1]);
    //adjust for the -1 of variant files for deletions
    auto it = std::find(tmp.nuc.begin(), tmp.nuc.end(), '-');
    if(it != tmp.nuc.end()){
      tmp.position = tmp.position+1;
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
      tmp.std_dev = std::stod(row_values[24]);
      if(row_values[26] != "NA"){
        tmp.amplicon_numbers = split_csv(row_values[26]);
      }
      if(row_values[25] != "NA"){
        tmp.freq_numbers = split_csv_double(row_values[25]);
      }
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

void set_insertion_flags(std::vector<variant> &variants){
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].depth_flag) continue;
    bool found = std::find(variants[i].nuc.begin(), variants[i].nuc.end(), '+') != variants[i].nuc.end();
    if(found){
      variants[i].include_clustering = false;
    }
  }
}

void set_deletion_flags(std::vector<variant> &variants, double lower_bound){
  std::vector<uint32_t> del_positions;
  std::vector<uint32_t> all_del_positions;

  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].depth_flag) continue;

    bool found = std::find(variants[i].nuc.begin(), variants[i].nuc.end(), '-') != variants[i].nuc.end();

    //here we divide by two because universal cluster could contain some noise of two var
    if(found && variants[i].gapped_freq > (lower_bound/(double)2 )){
      for(uint32_t j=1; j < variants[i].nuc.size()-1; j++){
        del_positions.push_back(variants[i].position+j);
      }
    }
    if(found && variants[i].gapped_freq > 0.001){
      for(uint32_t j=1; j < variants[i].nuc.size()-1; j++){
        all_del_positions.push_back(variants[i].position+j);
      }
    }
  }
  //set the include_clustering flag to false if this covers a deletion position
  for(uint32_t i=0; i < variants.size(); i++){
    bool found = std::find(del_positions.begin(), del_positions.end(), variants[i].position) != del_positions.end();
    size_t count = std::count(all_del_positions.begin(), all_del_positions.end(), variants[i].position);
    //if the positions is covered by two deletions we exclude it from clustering
    if(found || count > 1) {
      variants[i].include_clustering = false;
    } else {
      variants[i].include_clustering = true;
    }
  }
}


arma::mat form_dataset(std::vector<variant> base_variants, std::vector<variant> &variants, double lower_bound, double upper_bound){
  uint32_t useful_var=0;
  std::vector<uint32_t> count_pos;
  set_freq_range_flags(base_variants, lower_bound, upper_bound, true);
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked && base_variants[i].include_clustering){
      useful_var++;
      variants.push_back(base_variants[i]);
      count_pos.push_back(base_variants[i].position);
    }
  }

  arma::mat data(1, useful_var, arma::fill::zeros);

  //(rows, cols) where each columns is a sample
  for(uint32_t i = 0; i < variants.size(); i++){
    double tmp = static_cast<double>(variants[i].gapped_freq);
    data.col(i) = tmp;
  }
  return(data);

}

std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, std::vector<double> &solution, std::vector<double> &means, std::vector<double> &std_devs, std::string ref, double default_threshold, double &error_rate){
  if(ref.empty()){
    std::cerr << "Please provide a reference sequence." << std::endl;
    exit(1);
  }
  uint32_t n=8;
  uint32_t round_val = 4;
  bool development_mode=true;
  std::vector<variant> base_variants;
  parse_internal_variants(prefix, base_variants, min_depth, round_val, min_qual);
  set_deletion_flags(base_variants, 0.001);
  cluster_error(base_variants, min_qual, min_depth, error_rate);
  double lower_bound = 1-error_rate+0.0001;
  double upper_bound = error_rate-0.0001;
  set_freq_range_flags(base_variants, lower_bound, upper_bound, true);
  set_deletion_flags(base_variants, lower_bound);
  set_insertion_flags(base_variants);

  std::cerr << "lower bound " << lower_bound <<  " upper bound " << upper_bound << std::endl;

  uint32_t useful_var=0;
  std::vector<variant> variants;
  std::vector<uint32_t> count_pos;
  std::vector<uint32_t> second_dimension;
 
  set_freq_range_flags(base_variants, lower_bound, upper_bound, true);
  std::map<uint32_t, uint32_t> position_counts;

  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked && base_variants[i].include_clustering){
      useful_var++;
      if(position_counts.find(base_variants[i].position) != position_counts.end()){
        position_counts[base_variants[i].position] = 0;
      } else {
        position_counts[base_variants[i].position] += 1;
      }
      second_dimension.push_back(position_counts[base_variants[i].position]);
      variants.push_back(base_variants[i]);
      count_pos.push_back(base_variants[i].position);
      //std::cerr << "position " << base_variants[i].position << " nuc " << base_variants[i].nuc << " depth " << base_variants[i].depth << " gapped freq " << base_variants[i].gapped_freq <<  " include " << base_variants[i].include_clustering << std::endl;
    }
  }
  //std::cerr << "useful var " << useful_var << std::endl;

  if(useful_var < 1){
    std::ofstream file;
    if(development_mode){
      //write means to string
      file.open(output_prefix + ".txt", std::ios::trunc);
      std::string means_string = "[";
      means_string += "0.99";
      means_string += "]";
      file << "means\n";
      file << means_string << "\n";
      file.close();

      std::string solution_string = "[0.99]";
      std::string solution_filename = output_prefix + "_solution.txt";
      std::ofstream file_sol(solution_filename);
      file_sol << "means\n";
      file_sol << solution_string << "\n";
      file_sol.close();

    }
    call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
    return(variants);
  }

  uint32_t lower_n = find_max_frequency_count(count_pos);

  //initialize armadillo dataset and populate with frequency data
  arma::mat data(1, useful_var, arma::fill::zeros);

  //(rows, cols) where each columns is a sample
  for(uint32_t i = 0; i < variants.size(); i++){
    double tmp = static_cast<double>(variants[i].gapped_freq);
    data.col(i) = tmp;
  }

  arma::mat second_dim_data(2, useful_var, arma::fill::zeros);

  //(rows, cols) where each columns is a sample
  for(uint32_t i = 0; i < variants.size(); i++){
    double tmp = static_cast<double>(variants[i].gapped_freq);
    double perturb = (double)second_dimension[i];
    second_dim_data(0,i) = tmp;
    second_dim_data(1, i) = perturb;
  }

  //initialize things prior to clustering
  uint32_t counter = 1;
  uint32_t optimal_n = 0;
  gaussian_mixture_model retrained;

  std::vector<double> aics;
  std::vector<double> ns;
  std::vector<double> exclude_ns;
   while(counter <= n){
    //must have at least one point per cluster
    if(((double)useful_var < (double)counter)){
      break;
    }
    std::cerr << "\nn: " << counter << std::endl;
    retrained.means.clear();
    retrained.hefts.clear();
    retrained.prob_matrix.clear();
    bool clustering_failed = false;
    retrained = retrain_model(counter, data, second_dim_data, variants, lower_n, 0.001, clustering_failed);
    if(clustering_failed){
      counter++;
      continue;
    }

    std::vector<std::vector<double>> clusters = retrained.clusters;
    bool empty_cluster = false;

    std::vector<double> mads;
    for(auto data : clusters){
      double mean = calculate_mean(data);
      double mad = calculate_mad(data, mean);
      double std = calculate_standard_deviation(data);
      mads.push_back(mad);
      if(data.size() < 1){
        empty_cluster = true;
        std::cerr << "empty cluster " << counter << std::endl;
        if(counter == 2){
          optimal_n = 2;
          empty_cluster = false;
          break;
        }
        continue;
      }
      std::cerr << "mean " << mean << " mad " << mad << " std " << std << " cluster size " << data.size() <<  std::endl;
    }
    if(empty_cluster) {
      break;
    }
    //if the mean average deviation is low and the clusters are set to two, use this
    double threshold = 0.05;
    bool all_below = std::all_of(mads.begin(), mads.end(), [&](double v){ return v < threshold; });
    if(counter == 2 && all_below){
      optimal_n = counter;
      break;
    }
    threshold = 0.10;
    bool any_above = std::any_of(mads.begin(), mads.end(), [&](double v){ return v > threshold; });
    ns.push_back((double)counter);
    aics.push_back(retrained.bic);
    if(any_above){
      exclude_ns.push_back((double)counter);
    }
    counter++;
  }

  for(uint32_t i=0; i < aics.size(); i++){
    std::cerr << i+2 << " " << aics[i] << std::endl;
  }
  if(optimal_n == 0){
    std::cerr << "elbow gmm" << std::endl;
    optimal_n = elbow_method(aics, ns, exclude_ns);
  }
  std::cerr << "optimal n " << optimal_n << std::endl;


   if(optimal_n != retrained.means.size()){
    retrained.means.clear();
    retrained.hefts.clear();
    retrained.prob_matrix.clear();
    bool clustering_failed = false;
    retrained = retrain_model(optimal_n, data, second_dim_data, variants, lower_n, 0.001, clustering_failed);
  }
  for(auto cluster : retrained.clusters){
    double mean = calculate_mean(cluster);
    means.push_back(mean);
  }
  retrained.means = means;

  //TODO this can be put in function
  std::ofstream file;
  if(development_mode){
    //write means to string
    file.open(output_prefix + ".txt", std::ios::trunc);
    std::string means_string = "[";
    for(uint32_t j=0; j < retrained.means.size(); j++){
      if(j != 0){
        means_string += ",";
      }
      means_string += std::to_string(retrained.means[j]);
    }
    means_string += "]";
    file << "means\n";
    file << means_string << "\n";
    file.close();
  }

  assign_all_variants(variants, base_variants, retrained, lower_bound, upper_bound);
  add_noise_variants(variants, base_variants);
  if(retrained.n == 1){
    solution = retrained.means;
  }
  calculate_cluster_deviations(retrained);
  std_devs = retrained.cluster_std_devs;
  solve_clusters(variants, retrained, lower_bound, solution, output_prefix, default_threshold, min_depth);
  return(variants);
}
