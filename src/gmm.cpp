#include "gmm_1d.h"
#include "solve_clustering.h"
#include "genomic_position.h"
#include "gmm.h"
#include "saga.h"
#include "call_consensus_clustering.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

//weighted stdev of a variant's per-amplicon frequencies above which the position is
//masked. Defaults to 2.0, which is unreachable for frequencies bounded in [0,1], so
//amplicon masking is off unless a threshold is passed
const double DEFAULT_AMPLICON_STDEV = 2.0;

static double weighted_standard_deviation(const std::vector<double> &values, const std::vector<uint32_t> &weights){
  double weighted_sum = 0.0, total_weight = 0.0;
  for(uint32_t i=0; i < values.size(); i++){
    weighted_sum += values[i] * weights[i];
    total_weight += weights[i];
  }
  if(total_weight == 0) return 0.0;
  double mean = weighted_sum / total_weight;

  double variance = 0.0;
  for(uint32_t i=0; i < values.size(); i++){
    variance += weights[i] * std::pow(values[i] - mean, 2);
  }
  variance /= total_weight;
  return std::sqrt(variance);
}

//a variant is amplicon masked when an amplicon covering it fluctuates at some other
//position. Runs over every variant at once, so the result does not depend on the
//order positions were written in
void propagate_amplicon_masking(std::vector<variant> &variants){
  //count rather than flag, so a variant's own flux can be subtracted back out below
  std::unordered_map<std::string, uint32_t> flagged_counts;
  for(const auto &var : variants){
    if(!var.position_masked) continue;
    for(const auto &id : var.amplicon_ids){
      flagged_counts[id]++;
    }
  }

  for(auto &var : variants){
    uint32_t self = var.position_masked ? 1 : 0;
    var.amplicon_masked = false;
    for(const auto &id : var.amplicon_ids){
      auto it = flagged_counts.find(id);
      //a variant's own flux is reported by position_masked, not here
      if(it != flagged_counts.end() && it->second > self){
        var.amplicon_masked = true;
        break;
      }
    }
  }
}

void flag_amplicon_variation(std::vector<variant> &variants, double amplicon_stdev){
  for(auto &var : variants){
    if(var.freq_numbers.size() < 2 || var.freq_numbers.size() != var.amplicon_depths.size()){
      var.position_masked = false;
      continue;
    }
    var.position_masked = weighted_standard_deviation(var.freq_numbers, var.amplicon_depths) > amplicon_stdev;
  }
  propagate_amplicon_masking(variants);
}

void amplicon_specific_cluster_assignment(std::vector<variant> &variants, const gmm_1d &model, const std::vector<int> &component_indices){
  if(component_indices.empty()) return;
  for(uint32_t i=0; i < variants.size(); i++){
    variants[i].freq_assignments.clear();
    if(variants[i].freq_numbers.size() < 2) continue;
    if(!variants[i].position_masked) continue;

    //predict for all frequencies of this variant on each amplicon
    std::vector<std::vector<double>> proba = model.predict_proba(variants[i].freq_numbers);
    //per frequency assignment, restricted to the components the model kept so that
    //a discarded component cannot split two frequencies of the same cluster
    for(auto cluster_probs : proba){
      uint32_t best = component_indices[0];
      for(int ci : component_indices){
        if(cluster_probs[ci] > cluster_probs[best]) best = ci;
      }
      variants[i].freq_assignments.push_back(best);
    }
  }
}

void rewrite_position_masking(std::vector<variant> &variants){
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].freq_numbers.size() < 2) continue;
    if(!variants[i].position_masked) continue;
    if(variants[i].freq_assignments.empty()) continue;
    //per-amplicon frequencies that all land in one cluster are within-cluster noise,
    //not flux, so clear the stdev flag
    bool all_equal = std::all_of(variants[i].freq_assignments.begin(), variants[i].freq_assignments.end(), [&](uint32_t v) {return v == variants[i].freq_assignments[0];});
    if(all_equal){
      variants[i].position_masked = false;
    }
  }
  propagate_amplicon_masking(variants);
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

std::vector<std::string> split_csv_string(const std::string& input) {
    std::vector<std::string> result;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ',')) {
      result.push_back(token);
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

  //track which ref alleles we've already added
  while (std::getline(infile, line)) {
    std::vector<std::string> row_values;
    std::stringstream row_ss(line);
    std::string value;
    while (std::getline(row_ss, value, '\t')) {
      row_values.push_back(value);
    }
    //trailing empty fields are dropped by getline, so a named column may be absent
    //from the row even when the header declares it
    auto field = [&](const std::string &name) -> std::string {
      auto it = col_index.find(name);
      if(it == col_index.end() || (uint32_t)it->second >= row_values.size()) return "";
      return row_values[it->second];
    };
    variant tmp;
    tmp.nuc = field("ALT");
    tmp.position = std::stoi(field("POS"));
    //adjust for the -1 of variant files for deletions
    auto it = std::find(tmp.nuc.begin(), tmp.nuc.end(), '-');
    if(it != tmp.nuc.end()){
      tmp.position = tmp.position+1;
    }
    tmp.depth = std::stoi(field("ALT_DP"));
    tmp.total_depth = std::stoi(field("TOTAL_DP"));
    tmp.freq = std::round(std::stof(field("ALT_FREQ")) * multiplier) / multiplier;
    tmp.qual = std::stod(field("ALT_QUAL"));

    if(col_index.count("GAPPED_FREQ")){
      tmp.gapped_freq = round(std::stod(field("GAPPED_FREQ")) * multiplier) / multiplier;
      tmp.gapped_depth = std::stoi(field("GAPPED_DEPTH"));
      if(!is_empty_field(field("AMP_NUMBERS"))){
        tmp.amplicon_ids = split_csv_string(field("AMP_NUMBERS"));
      }
      if(!is_empty_field(field("AMP_FREQ"))){
        tmp.freq_numbers = split_csv_double(field("AMP_FREQ"));
      }
      if(!is_empty_field(field("AMP_DEPTH"))){
        tmp.amplicon_depths = split_csv(field("AMP_DEPTH"));
      }
      tmp.version_1_var = false;
    } else {
      tmp.gapped_freq = 0.0;
      tmp.gapped_depth = 0;
      tmp.version_1_var = true;
    }
    if(tmp.gapped_depth < depth_cutoff){
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
      //variants[i].include_clustering = false;
    }
  }
}

void set_deletion_flags(std::vector<variant> &variants, double lower_bound, double invariant_lower_bound){
  // Collect all valid deletions with their covered ranges.
  struct del_info {
    uint32_t idx;
    uint32_t start;
    uint32_t end;
    double freq;
  };

  std::vector<del_info> dels;
  for (uint32_t i = 0; i < variants.size(); i++) {
    if (variants[i].depth_flag) continue;
    if (variants[i].nuc.find('-') == std::string::npos) continue;
    std::string nuc = variants[i].nuc;
    nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
    if (nuc.empty()) continue;
    uint32_t end = variants[i].position + (uint32_t)nuc.size() - 1;
    dels.push_back({i, variants[i].position, end, variants[i].gapped_freq});
  }

  // Sort dominant-first; greedily accept non-overlapping deletions and mark
  // any deletion that overlaps an already-accepted one as overlapped_deletion.
  std::sort(dels.begin(), dels.end(), [](const del_info& a, const del_info& b) {
    return a.freq > b.freq;
  });

  std::vector<std::pair<uint32_t, uint32_t>> accepted;
  for (const auto& d : dels) {
    bool overlaps = false;
    for (const auto& [s, e] : accepted) {
      if (d.start <= e && d.end >= s) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) {
      variants[d.idx].overlapped_deletion = true;
    } else {
      accepted.push_back({d.start, d.end});
    }
  }

  // Accumulate gap signal from sub-threshold deletions so that variants whose
  // reduced gapped_freq is fully explained by noise-level deletions are treated
  // as invariant rather than assigned to a single cluster.
  std::unordered_map<uint32_t, double> sub_threshold_del_freq;
  for (uint32_t i = 0; i < variants.size(); i++) {
    if (variants[i].depth_flag) continue;
    if (variants[i].nuc.find('-') == std::string::npos) continue;
    if (!variants[i].outside_freq_range) continue;
    if (variants[i].gapped_freq > invariant_lower_bound) continue;
    std::string nuc = variants[i].nuc;
    nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
    for (uint32_t j = 1; j < nuc.size(); j++) {
      sub_threshold_del_freq[variants[i].position + j] += variants[i].gapped_freq;
    }
  }

  double invariant_upper_bound = 1.0 - invariant_lower_bound;
  for (uint32_t i = 0; i < variants.size(); i++) {
    if (variants[i].outside_freq_range) continue;
    auto it = sub_threshold_del_freq.find(variants[i].position);
    if (it == sub_threshold_del_freq.end()) continue;
    if (variants[i].gapped_freq + it->second >= invariant_upper_bound) {
      variants[i].half_normal_upper = true;
      variants[i].outside_freq_range = true;
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
                              double default_threshold, \
                              uint32_t n, double invariant_threshold, double covariance_prior, double mean_precision_prior, double min_cluster_fraction, double amplicon_stdev){
  uint32_t round_val = 4;
  std::vector<variant> base_variants;
  parse_internal_variants(prefix, base_variants, min_depth, round_val, min_qual, invariant_threshold);
  flag_amplicon_variation(base_variants, amplicon_stdev);
  set_deletion_flags(base_variants, 0.001, 1.0 - invariant_threshold);

  std::vector<variant> model_variants;
  std::vector<double> model_freqs;
  std::vector<double> all_freqs;
  for(uint32_t i=0; i < base_variants.size(); i++){
    all_freqs.push_back(base_variants[i].gapped_freq);
    if(!base_variants[i].depth_flag && !base_variants[i].qual_flag && !base_variants[i].outside_freq_range && !base_variants[i].overlapped_deletion){
      model_variants.push_back(base_variants[i]);
      model_freqs.push_back(base_variants[i].gapped_freq);
      //std::cerr << base_variants[i].nuc << " " << base_variants[i].position << " " << base_variants[i].gapped_freq << "\n";
    }
  }
  std::cerr << "Number of frequencies used for clustering: " << model_freqs.size() << "\n";
  
  //handle the case of no variants less than the universal cluster
  if(model_variants.size() <= 1){
    variant_assigner(std::vector<double>(1, 1.0), std::vector<double>(1, 1.0), 2.0).assign(base_variants);
    call_majority_consensus(base_variants, output_prefix, default_threshold);
    write_single_cluster_output(output_prefix);
    base_variants.clear();
    return(base_variants);
  }
  gmm_1d model(n, 42);  // seed matches bootstrap replicate
  model.set_use_half_normal_for_noise(true, invariant_threshold);

  model.set_covariance_prior(covariance_prior);
  model.set_mean_precision_prior(mean_precision_prior);

  model.fit(model_freqs);
  model.set_min_cluster_fraction(min_cluster_fraction);

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
 
  std::cerr << "all freqs " << all_freqs.size() << " base variants " << base_variants.size() << "\n";
  //take the model labels and assign them to the variants
  for(uint32_t i=0; i < all_labels.size(); i++){
    base_variants[i].cluster_assigned = all_labels[i];
    if(!base_variants[i].half_normal_upper && !base_variants[i].half_normal_lower) {
      base_variants[i].probabilities = proba[i];
    }

    if(base_variants[i].gapped_freq > invariant_threshold){
      base_variants[i].half_normal_upper = true;
      base_variants[i].half_normal_lower = false;
    } else if(base_variants[i].gapped_freq < 1-invariant_threshold){
      base_variants[i].half_normal_lower = true;
      base_variants[i].half_normal_upper = false;
    } else if(all_labels[i] == n-1){
      base_variants[i].half_normal_upper = true;
      base_variants[i].half_normal_lower = false;
    } else if(all_labels[i] == n-2){
      base_variants[i].half_normal_lower = true;
      base_variants[i].half_normal_upper = false;
    }

  }
  
  subset_sum_solver solver(eff_means, subset_sum_solver::UNIT_SUM_ERROR, invariant_threshold);
  bool solved = solver.solve();
  //the boundary rescue may have refined a mean, and downstream matching against the
  //solution is exact equality, so adopt the refined means
  eff_means = solver.refined_means();
  std::vector<std::vector<double>> solution_sets = solver.get_solution_sets();

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
      variant_assigner(std::vector<double>(1, 1.0), std::vector<double>(1, 1.0), 2.0).assign(base_variants);
      call_majority_consensus(base_variants, output_prefix, default_threshold);
      base_variants.clear();
    } else{
      variant_assigner::overwrite_cluster_assigned(base_variants, eff_means, model_means);
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
      variant_assigner(solution_sets[0], eff_means, 2.0).assign(base_variants);

      //predict clusters for amplicon specific frequencies
      amplicon_specific_cluster_assignment(base_variants, model, component_indices);
      //write the amplicon flags based on cluster agreement
      rewrite_position_masking(base_variants);

      solution = solution_sets[0];
    }
  } else {
    variant_assigner(std::vector<double>(1, 1.0), std::vector<double>(1, 1.0), 2.0).assign(base_variants);
    call_majority_consensus(base_variants, output_prefix, default_threshold);
    base_variants.clear();
  }  

  means = eff_means;
  return(base_variants);
}