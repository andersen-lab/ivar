#include "gmm_1d.h"
#include "solve_clustering.h"
#include "genomic_position.h"
#include "gmm.h"
#include "saga.h"
#include "call_consensus_clustering.h"
#include "ref_seq.h"
#include <fstream>
#include <cassert>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

void flag_low_posterior_variants(std::vector<variant> &base_variants){
  for(auto &var : base_variants){
    if(var.probabilities.empty()) continue;
    uint32_t assigned_cluster = var.cluster_assigned;
    double assigned_prob = var.probabilities[assigned_cluster];

    double next_best = 0.0;
    for(uint32_t k = 0; k < var.probabilities.size(); k++){
      if(k == assigned_cluster) continue;
      if(var.probabilities[k] > next_best) next_best = var.probabilities[k];
    }

    if(next_best > 0.0 && (assigned_prob / next_best) < 2.0){
      var.vague_assignment = true;
    } else {
      var.vague_assignment = false;
    }
  }
}

void flag_position_conflicts(std::vector<variant> &variants) {
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> pos_cluster_count;

  for (const auto& v : variants) {
    if (v.half_normal_upper || v.half_normal_lower || v.depth_flag || v.qual_flag) continue;
    /*if(v.position == 29402){
      std::cerr << "position " << v.position << " " << v.gapped_freq << " " << v.cluster_assigned << std::endl;
    }*/
    pos_cluster_count[v.position][v.cluster_assigned]++;
  }

  std::unordered_set<uint32_t> conflicted;
  for (const auto& [pos, cluster_counts] : pos_cluster_count) {
    for (const auto& [cluster, count] : cluster_counts) {
      if (count > 1) {
        conflicted.insert(pos);
        break;
      }
    }
  }

  for (auto& v : variants) {
    if (conflicted.count(v.position)) {
      v.position_conflict = true;
    }
  }
}

void reset_variants_info(std::vector<variant> &variants){
  //reset cluster assignments and probabilities prior to rerunning another model
  for(auto &var : variants){
    var.cluster_assigned = -1;
    var.probabilities.clear();
  }
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

std::vector<double> loglikelihoods_to_posteriors(const std::vector<double>& loglikes){
    const size_t K = loglikes.size();
    if (K == 0) return {};
    // 1. Find max for numerical stability
    double m = *std::max_element(loglikes.begin(), loglikes.end());
    // 2. Exponentiate shifted log-likelihoods
    std::vector<double> exps(K);
    double sum_exp = 0.0;
    for (size_t k = 0; k < K; ++k) {
        exps[k] = std::exp(loglikes[k] - m);
        sum_exp += exps[k];
    }
    // 3. Normalize
    for (size_t k = 0; k < K; ++k) {
        exps[k] /= sum_exp;
    }
    return exps;
}

std::vector<uint32_t> compare_cluster_assignment(std::vector<std::vector<double>> prob_matrix, std::vector<uint32_t> assigned){
  double threshold = 2;
  std::vector<uint32_t> flagged_idx;

  for(uint32_t i=0; i < prob_matrix.size(); i++){
    std::vector<double> probs = loglikelihoods_to_posteriors(prob_matrix[i]);
    /*for(auto p : probs){
      std::cerr << p << " ";
    }
    std::cerr << "\n";*/
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

bool is_empty_field(const std::string& val) {
    if (val.empty()) return true;

    // check if only whitespace
    if (val.find_first_not_of(" \t\n\r") == std::string::npos)
        return true;

    for (unsigned char c : val) {
      if (!(std::isprint(c) || std::isspace(c))) {
        return true;
      }
    }

    // common missing-value markers
    if (val == "NA" || val == "null" || val == "None")
        return true;

    return false;
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

void parse_internal_variants(std::string filename, 
  std::vector<variant> &variants, 
  uint32_t depth_cutoff, 
  uint32_t round_val, 
  uint8_t quality_threshold, 
  double invariant_threshold){

  std::ifstream infile(filename);
  std::string line;

  // --- Read header ---
  std::getline(infile, line);
  std::stringstream ss(line);
  std::string col;
  std::vector<std::string> headers;
  std::unordered_map<std::string, int> col_index;

  int idx = 0;
  while (std::getline(ss, col, '\t')) {
    headers.push_back(col);
    col_index[col] = idx++;
  }

  uint32_t count = 0;
  double multiplier = pow(10, round_val);
  double compare_quality = static_cast<double>(quality_threshold);

  auto to_bool = [](const std::string& s) -> bool {return s == "TRUE" || s == "true" || s == "True" || s == "1";};
  //track which ref alleles we've already added
  while (std::getline(infile, line)) {
    std::vector<std::string> row_values;
    std::stringstream row_ss(line);
    std::string value;
    while (std::getline(row_ss, value, '\t')) {
      row_values.push_back(value);
    }
    variant tmp;
    tmp.nuc = row_values[col_index["ALT"]];
    tmp.position = std::stoi(row_values[col_index["POS"]]);
    //adjust for the -1 of variant files for deletions
    auto it = std::find(tmp.nuc.begin(), tmp.nuc.end(), '-');
    if(it != tmp.nuc.end()){
      tmp.position = tmp.position+1;
    }
    tmp.depth = std::stoi(row_values[col_index["ALT_DP"]]);
    tmp.total_depth = std::stoi(row_values[col_index["TOTAL_DP"]]);
    tmp.freq = std::round(std::stof(row_values[col_index["ALT_FREQ"]]) * multiplier) / multiplier;
    tmp.qual = std::stod(row_values[col_index["ALT_QUAL"]]);

    if(row_values.size() > 20){
      tmp.gapped_freq = round(std::stod(row_values[col_index["GAPPED_FREQ"]]) * multiplier) / multiplier;
      tmp.gapped_depth = std::stoi(row_values[col_index["GAPPED_DEPTH"]]);
      tmp.amplicon_flux = to_bool(row_values[col_index["FLAGGED_POS"]]);
      tmp.amplicon_masked = to_bool(row_values[col_index["AMP_MASKED"]]);
      if(!(is_empty_field(row_values[col_index["STD_DEV"]]))){
        tmp.std_dev = std::stod(row_values[col_index["STD_DEV"]]);
      } else {
        tmp.std_dev = 0;
      }
      if(row_values.size() > 26 && !(is_empty_field(row_values[col_index["AMP_NUMBERS"]]))){ 
        tmp.amplicon_numbers = split_csv(row_values[col_index["AMP_NUMBERS"]]);
      } else {
        tmp.amplicon_numbers = {};
      }
      if(row_values.size() > 25 && !(is_empty_field(row_values[col_index["AMP_FREQ"]]))){
        tmp.freq_numbers = split_csv_double(row_values[col_index["AMP_FREQ"]]);
      } else {
        tmp.freq_numbers = {};
      }
      tmp.version_1_var = false;
    } else {
      tmp.gapped_freq = 0.0;
      tmp.amplicon_flux = false;
      tmp.amplicon_masked = false;
      tmp.primer_masked = false;
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
    if(tmp.gapped_freq >= invariant_threshold || tmp.gapped_freq <= (1-invariant_threshold)){
      tmp.outside_freq_range = true;
    } else {
      tmp.outside_freq_range = false;
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

void write_single_cluster_output(std::string output_prefix){
  std::ofstream out(output_prefix + "_gmm_1d_results.txt");
  out << "Components\tDistinct_Components\tMeans\tVariances\tWeights\tEffective_Means\tEffective_Variances\tEffective_Weights\tSolution_Sets\n";
  out << "1\t1\t[]\t[]\t[]\t[0,1]\t[0,0]\t[0,0]\t[[1]]\n";
  out.close();
}

std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, \
                              std::vector<double> &solution, std::vector<double> &means, \
                              std::string ref, double default_threshold, \
                              uint32_t n, double invariant_threshold, double covariance_prior, double mean_precision_prior){
  if(ref.empty()){
    std::cerr << "Please provide a reference sequence." << std::endl;
    exit(1);
  }

  uint32_t round_val = 4;
  std::vector<variant> base_variants;
  parse_internal_variants(prefix, base_variants, min_depth, round_val, min_qual, invariant_threshold);
  set_deletion_flags(base_variants, 0.001);

  std::vector<variant> model_variants;
  std::vector<double> model_freqs;
  std::vector<double> all_freqs;
  for(uint32_t i=0; i < base_variants.size(); i++){
    all_freqs.push_back(base_variants[i].gapped_freq);
    if(!base_variants[i].depth_flag && !base_variants[i].qual_flag && !base_variants[i].outside_freq_range){
      model_variants.push_back(base_variants[i]);
      model_freqs.push_back(base_variants[i].gapped_freq);
      //std::cerr << base_variants[i].position << " " << base_variants[i].gapped_freq << "\n";
    }
  }
  std::cerr << "Number of frequencies used for clustering: " << model_freqs.size() << "\n";
  
  //handle the case of no variants less than the universal cluster
  if(model_variants.size() <= 1){
    call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
    write_single_cluster_output(output_prefix);
    base_variants.clear();
    return(base_variants);
  }
  gmm_1d model(n, 42);  // seed matches bootstrap replicate
  model.set_use_half_normal_for_noise(true, invariant_threshold);

  // spike in priors
  //model.set_mean_precision_prior(0.5);

  //simulated priors
  //model.set_covariance_prior(1e-3);
  //model.set_mean_precision_prior(1e-2);

  model.set_covariance_prior(covariance_prior);
  model.set_mean_precision_prior(mean_precision_prior);

  model.fit(model_freqs);
  model.set_min_cluster_fraction(0.10);

  std::vector<int> labels = model.predict(model_freqs);

  std::cerr << "All means: ";
  for(auto x: model.get_means()){
    std::cerr << x << " ";
  }
  std::cerr << "\n";

  std::vector<int> component_indices = model.get_effective_components(labels);
  std::cerr << "VB effective components: " << component_indices.size() << "\n";
  std::cerr << "VB effective means: ";
    std::vector<double> eff_means = model.get_effective_means(component_indices);
    for(int i = 0; i < eff_means.size(); i++) {
      std::cerr << eff_means[i] << " ";
    }
    std::cerr << "\n";

    std::cerr << "VB effective vars: ";
    std::vector<double> eff_vars = model.get_effective_vars(component_indices);
    for(int i = 0; i < eff_vars.size(); i++) {
      std::cerr << eff_vars[i] << " ";
    }
    std::cerr << "\n";

    std::cerr << "VB effective weights: ";
    std::vector<double> eff_weights = model.get_effective_weights(component_indices);
    for(int i = 0; i < eff_weights.size(); i++) {
      std::cerr << eff_weights[i] << " ";
    }
    std::cerr << "\n";

    std::cerr << model.get_invariant_components(labels).size() << " invariant components\n";

    std::cerr << "ELBO history: ";
    for(int i = 0; i < model.get_elbo_history().size(); i++) {
      std::cerr << model.get_elbo_history()[i] << ", ";
    }
    std::cerr << "\n";

  //predict on all the frequencies - disable min_cluster_fraction here because all_freqs
  //includes tens of thousands of reference alleles, making the 10% threshold nonsensical
  //for the small signal clusters
  model.set_min_cluster_fraction(0.0);
  std::vector<int> all_labels = model.predict(all_freqs);
  //gets the posterior probability per variant
  std::vector<std::vector<double>> proba = model.predict_proba(all_freqs);
  
  std::vector<double> model_means = model.get_means();
  std::vector<double> model_vars = model.get_variances();
  std::vector<double> model_weights = model.get_weights();

  //take the model labels and assign them to the variants
  for(uint32_t i=0; i < all_labels.size(); i++){
    if(all_labels[i] == n-1 || base_variants[i].gapped_freq > invariant_threshold){
      base_variants[i].half_normal_upper = true;
      base_variants[i].half_normal_lower = false;
    } else if(all_labels[i] == n-2 || base_variants[i].gapped_freq < 1-invariant_threshold){
      base_variants[i].half_normal_lower = true;
      base_variants[i].half_normal_upper = false;
    }
    base_variants[i].cluster_assigned = all_labels[i];
    if(!base_variants[i].half_normal_upper && !base_variants[i].half_normal_lower) {
      base_variants[i].probabilities = proba[i];
    }
  }
  
  std::vector<std::vector<double>> solution_sets;
  bool solved = subset_sum(eff_means, solution_sets, 0.10);
 

  //output the clustering information
  std::ofstream out(output_prefix + "_gmm_1d_results.txt");
  if (!out.is_open()) {
    std::cerr << "Failed to open output file: " << output_prefix << std::endl;
  }
  if (!out) {
    std::cerr << "Stream bad immediately after open: " << output_prefix << std::endl;
  }
  out << "Components\tDistinct_Components\tMeans\tVariances\tWeights\tEffective_Means\tEffective_Variances\tEffective_Weights\tSolution_Sets\n";
  out << std::to_string(n) << "\t";
  out << component_indices.size() << "\t";
  out << "[";
    for (size_t j = 0; j < model_means.size(); j++) {
      out << model_means[j];
      if (j < means.size() - 1) out << ",";
    }
  out << "]";
  out << "\t";
  out << "[";
  for (size_t j = 0; j < model_vars.size(); j++) {
    out << model_vars[j];
    if (j < model_vars.size() - 1) out << ",";
  }
  out << "]";
  out << "\t";

  out << "[";
  for (size_t j = 0; j < model_weights.size(); j++) {
    out << model_weights[j];
    if (j < model_weights.size() - 1) out << ",";
  }
  out << "]";
  out << "\t";
  out << "[";
  for (size_t j = 0; j < eff_means.size(); j++) {
    out << eff_means[j];
    if (j < eff_means.size() - 1) out << ",";
  }
  out << "]";
  out << "\t";
  out << "[";
  for (size_t j = 0; j < eff_vars.size(); j++) {
    out << eff_vars[j];
    if (j < eff_vars.size() - 1) out << ",";
  }
  out << "]";
  out << "\t";

  out << "[";
  for (size_t j = 0; j < eff_weights.size(); j++) {
    out << eff_weights[j];
    if (j < eff_weights.size() - 1) out << ",";
  }
  out << "]";
  out << "\t";
  out << "[";

  for(uint32_t t=0; t < solution_sets.size(); t++){
    auto sol = solution_sets[t];
    if(t != 0) out << ",";
    out << "[";
    for(uint32_t s=0; s < sol.size(); s++){
      if(s != 0) out << ",";
      out << std::to_string(sol[s]);
    }
    out << "]";
  }
  out << "]";
  out << "\n";

  if(solved){
    if(solution_sets.size() > 1){
      call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
      base_variants.clear();
    } else{
      overwrite_cluster_assigned(base_variants, eff_means, model_means);

      //recalculate probabilities based on the new cluster assignments and only the effective means
      for(auto &v : base_variants){
        if(v.half_normal_upper || v.half_normal_lower || v.probabilities.empty()) continue;
        std::vector<double> eff_proba;
        double sum = 0.0;
        for(int ci : component_indices){
          eff_proba.push_back(v.probabilities[ci]);
          sum += v.probabilities[ci];
        }
        if(sum > 0.0)
          for(auto &p : eff_proba) p /= sum;
        v.probabilities = eff_proba;
      }
      //flag low posterior variants based on the new probabilities
      flag_low_posterior_variants(base_variants);
      assign_variants_solution(solution_sets[0], base_variants, eff_means);
      flag_position_conflicts(base_variants);
      solution = solution_sets[0];
    }
  } else {
    call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
    base_variants.clear();
  }  

  means = eff_means;
  return(base_variants);
}
