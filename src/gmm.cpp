#include "./include/armadillo"
#include "solve_clustering.h"
#include "gmm.h"
#include "saga.h"
#include "call_consensus_clustering.h"
#include "estimate_error.h"
#include "ref_seq.h"
#include <fstream>
#include <cassert>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>

void reset_variants_info(std::vector<variant> &variants){
  //reset cluster assignments and probabilities prior to rerunning another model
  for(auto &var : variants){
    var.cluster_assigned = -1;
    var.probabilities.clear();
  }
}

arma::mat subsample_with_replacement(
    const arma::mat& data,
    std::size_t n_subsample,
    const std::vector<uint32_t>& position,
    std::vector<variant> &subsampled_variants,
    std::vector<variant> variants) {
    
    // Build mapping: group_id -> indices of columns belonging to it
    std::unordered_map<uint32_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < position.size(); ++i) {
        groups[position[i]].push_back(i);
    }

    // Collect unique group IDs into a vector (for random selection)
    std::vector<uint32_t> group_ids;
    group_ids.reserve(groups.size());
    for (const auto& kv : groups) {
        group_ids.push_back(kv.first);
    }

    // Random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> dist(0, group_ids.size() - 1);

    std::vector<arma::uword> chosen_cols;
    chosen_cols.reserve(n_subsample); // approximate

    uint32_t pos = 0;
    for (std::size_t i = 0; i <= n_subsample; ++i) {
      uint32_t group = group_ids[dist(gen)];
      const auto& cols = groups[group];
      chosen_cols.insert(chosen_cols.end(), cols.begin(), cols.end());
      for(auto c : cols){
        auto &var = variants[c];
        var.position = pos;
        subsampled_variants.push_back(var);
      }
      i += cols.size();
      i--;
      pos += 1;
    }

    // Build subsample matrix from chosen columns
    arma::mat subsample = data.cols(arma::uvec(chosen_cols));

    return subsample;
}

double calculate_BIC(double k, double logL, int N) {
    if (N == 0) return std::numeric_limits<double>::infinity();
    double bic = -2.0 * logL + k * std::log(N);
    return bic;
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

  //assignment by position
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
      //std::cerr << "HERE " << tmp_prob.size() << " " << assigned.size() << " " << pos << std::endl;
      clustering_failed = true;
      continue;
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
    if((!base_variants[i].outside_freq_range && (base_variants[i].amplicon_flux || base_variants[i].amplicon_masked || base_variants[i].imbalance || !base_variants[i].include_clustering)) || (base_variants[i].outside_freq_range && base_variants[i].gapped_freq >= lower_bound && base_variants[i].gapped_freq <= upper_bound)){
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
}

void add_noise_variants(std::vector<variant> &variants, std::vector<variant> base_variants){
  //lets add back in the 100% variants
  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].outside_freq_range){
      variants.push_back(base_variants[i]);
    }
  }
}

//test function
std::string vec_to_pylist(const std::vector<double>& v){
    std::ostringstream ss;
    ss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i + 1 < v.size()) ss << ", ";
    }
    ss << "]";
    return ss.str();
}

gaussian_mixture_model retrain_model(uint32_t n, 
                                    arma::mat data,
                                    std::vector<variant> &variants, 
                                    uint32_t lower_n, 
                                    double &var_floor, 
                                    bool &clustering_failed, 
                                    bool error_clustering){

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

  gaussian_mixture_model gmodel;
  gmodel.n = n;
  gmodel.lower_n = lower_n;
  arma::gmm_diag model;

  if(error_clustering){
    var_floor = 0.0001;
  } else {
    var_floor = var_floor;
  }

  bool status = model.learn(data, n, arma::eucl_dist, arma::static_spread, 10, 10, var_floor, false);
  if(!status){
    std::cerr << "GMM failed to converge" << std::endl;
    clustering_failed = true;
    return(gmodel);
  }
  std::vector<double> means;
  
  for(auto h: model.hefts){
    double heft = (double)h;
    gmodel.hefts.push_back(heft);
  }

  for(auto d : model.dcovs){
    gmodel.dcovs.push_back(std::sqrt((double)d));
  }

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
  for(uint32_t i=0; i < model.means.size(); i++){
    double m = (double)model.means[i];
    means.push_back(m);
  }

  gmodel.prob_matrix = prob_matrix;
  gmodel.model = model;
  gmodel.means = means;

  std::vector<std::vector<uint32_t>> possible_permutations;
  for (uint32_t i = 1; i <= gmodel.lower_n; ++i){
    perm_generator(gmodel.n, i, possible_permutations);
  }

  assign_clusters(variants, gmodel, clustering_failed, possible_permutations, unique_pos, pos_to_variant_indices);
  std::vector<std::vector<double>> clusters = form_clusters(n, variants);
  gmodel.clusters = clusters;
  means.clear();
  //set means
  arma::mat mean_fill2 (1, n, arma::fill::zeros);
  for(uint32_t i=0; i < clusters.size(); i++){
    double m = calculate_mean(clusters[i]);
    double factor = std::pow(10.0, 2);
    double rounded = std::round(m * factor) / factor;
    mean_fill2.col(i) = rounded;
    means.push_back(rounded);

  }
  model.set_means(mean_fill2);
  gmodel.means = means;
  double k = (2 * n) + (n-1);
  double bic = calculate_BIC(k, gmodel.model.sum_log_p(data), (int) data.n_cols);
  gmodel.bic = bic;
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

std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, \
                              std::vector<double> &solution, std::vector<double> &means, std::vector<double> &std_devs, \
                              std::string ref, double default_threshold, double &error_rate){
  if(ref.empty()){
    std::cerr << "Please provide a reference sequence." << std::endl;
    exit(1);
  }
  uint32_t n=10;
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
  std::cerr << "filename: " << output_prefix << std::endl;
  std::cerr << "lower bound " << lower_bound <<  " upper bound " << upper_bound << std::endl;

  uint32_t useful_var=0;
  std::vector<variant> variants;
  std::vector<uint32_t> count_pos; 
  set_freq_range_flags(base_variants, lower_bound, upper_bound, true);

  std::unordered_map<uint32_t, uint32_t> position_counts;
  // first pass: count how many times each position appears
  for (uint32_t i=0; i < base_variants.size(); i++) {
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked && base_variants[i].include_clustering){ 
      position_counts[base_variants[i].position]++;
    }
  }

  for(uint32_t i=0; i < base_variants.size(); i++){
    if (position_counts[base_variants[i].position] < 2) {
      base_variants[i].imbalance = true;
      continue;
    }
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked && base_variants[i].include_clustering){
      useful_var++;
      variants.push_back(base_variants[i]);
      count_pos.push_back(base_variants[i].position);
      std::cerr << "position " << base_variants[i].position << " nuc " << base_variants[i].nuc << " depth " << base_variants[i].depth << " gapped freq " << base_variants[i].gapped_freq <<  " include " << base_variants[i].include_clustering << std::endl;
    }
  }
  //handle the case of no variants less than the universal cluster
  if(useful_var < 1){
    std::ofstream file;
    if(development_mode){
      //write means to string
      file.open(output_prefix + ".txt", std::ios::trunc);
      std::string means_string = "[[";
      means_string += "0.99";
      means_string += "]]";
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
  std::vector<uint32_t> subsample_position;
  for(uint32_t i = 0; i < variants.size(); i++){
    double tmp = static_cast<double>(variants[i].gapped_freq);
    subsample_position.push_back(variants[i].position);
    data.col(i) = tmp;
  }
  std::cerr << "useful var " << useful_var << std::endl;

  std::unordered_map<uint32_t, uint32_t> model_counter; 

  bool empty_cluster = false;
  std::vector<std::vector<double>> solution_sets;

  uint32_t counter = 1;
  bool clustering_failed =false;
  std::vector<variant> subsampled_variants;
  uint32_t bootstrap_reps = 1000;
  uint32_t final_n=0;
  double var_floor;
  std::vector<double> all_var_floors = {0.005}; 

  std::vector<double> track_var_floors;
  std::vector<double> track_bics;
  std::vector<std::vector<double>> track_means;
  std::vector<std::vector<double>> track_std_devs;
  std::vector<uint32_t> track_ns;
  std::vector<uint32_t> track_bootstraps;

  for(uint32_t i =0; i < bootstrap_reps; i++){
    empty_cluster = false;
    subsampled_variants.clear();
    arma::mat subsample = subsample_with_replacement(data, data.size(), subsample_position, subsampled_variants, variants);
    counter = 1;
    std::vector<double> all_bics;
    for(uint32_t j=0; j < all_var_floors.size(); j++){
      counter = 1;
      while(counter <= n && !empty_cluster){ 
          clustering_failed = false;
          reset_variants_info(subsampled_variants);
          reset_variants_info(variants);
          var_floor = all_var_floors[j];
          gaussian_mixture_model retrained = retrain_model(counter, subsample, subsampled_variants, lower_n, var_floor, clustering_failed, false);
          if(clustering_failed){
            all_bics.push_back(std::numeric_limits<double>::max());
            counter++;
            continue;  
          }

          //calculate_cluster_deviations(retrained);
          track_ns.push_back(counter);
          track_var_floors.push_back(var_floor);
          track_means.push_back(retrained.means);
          track_std_devs.push_back(retrained.dcovs);
          track_bics.push_back(retrained.bic);
          track_bootstraps.push_back(i+1);

          std::cerr << "bootstrap rep " << i << " bic " << retrained.bic << " n " << counter << "  var floor " << var_floor << std::endl;
          double bic = retrained.bic;
          all_bics.push_back(bic);
          counter++;
        }
    }
    double lowest = std::numeric_limits<double>::max(); 
    uint32_t winner;
    for(uint32_t i=0; i < all_bics.size(); i++){
      if(all_bics[i] < lowest){
        lowest = all_bics[i];
        winner = i+1;
      }
      std::cerr << i+1 << " " << all_bics[i] << "\n";
    }
    model_counter[winner] += 1;
  }

  uint32_t highest = 0;
  for (const auto& [key, value] : model_counter) {
    if(value >= 900){
      final_n = key;
      highest = value;
    }
    std::cerr << key << " : " << value << '\n';
  }
  gaussian_mixture_model best_model; 
  std::cerr << "final n " << final_n << std::endl;
  if(final_n ==0){
    std::cerr << output_prefix << " no solution found" << std::endl;
    exit(1);
  }

  std::cerr << "new final n " << final_n << std::endl;


  clustering_failed = false;
  reset_variants_info(variants);
  best_model= retrain_model(final_n, data, variants, lower_n, var_floor, clustering_failed, false); 
  calculate_cluster_deviations(best_model);
  
  for(auto m : best_model.means){
    std::cerr << "mean " << m << std::endl;  
  }

  //n, means, output_prefix, bic, var_floor, widths 
  /*std::ofstream out("bootstrap_stats_0.005_0.90.tsv", std::ios::app);
  for(uint32_t i=0; i < track_bics.size(); i++){
    out << track_ns[i] << "\t"
        << vec_to_pylist(track_means[i]) << "\t"
        << output_prefix << "\t"
        << std::to_string(track_bics[i]) << "\t"
        << std::to_string(track_var_floors[i]) << "\t"
        << vec_to_pylist(track_std_devs[i]) << "\t"
        << std::to_string(track_bootstraps[i]) << "\n";
  }
  out.close();*/
  
  std::ofstream file;
  if(development_mode){
    //write means to string
    file.open(output_prefix + ".txt", std::ios::trunc);
    std::string means_string = "[";
    for(uint32_t j=0; j < best_model.means.size(); j++){
      if(j != 0){
        means_string += ",";
      }
      means_string += std::to_string(best_model.means[j]);
    }
    means_string += "]";
    file << "means\n";
    file << means_string << "\n";
    file.close();
  }
  assign_all_variants(variants, base_variants, best_model, lower_bound, upper_bound);
  add_noise_variants(variants, base_variants);
  solve_clusters(variants, best_model, lower_bound, solution, output_prefix, default_threshold, min_depth);
  means = best_model.means;
  return(variants);
}
