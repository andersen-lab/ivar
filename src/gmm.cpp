#include "solve_clustering.h"
#include "gmm.h"
#include "gmm_1d.h"
#include "saga.h"
#include "call_consensus_clustering.h"
#include "estimate_error.h"
#include "ref_seq.h"
#include <numeric>
#include <fstream>
#include <cassert>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>

void compare_component_assigments(std::vector<variant> &base_variants, const double threshold=1.25){
  for(uint32_t i=0; i < base_variants.size(); i++){
    variant &var = base_variants[i];
    std::vector<double> tmp = var.marginal_posterior_probabilities;
    if(tmp.size() == 0) continue;

    std::sort(tmp.begin(), tmp.end(), std::greater<double>());
    double ratio = tmp[0] / tmp[1];
    if(ratio < threshold){
      var.vague_component_assignment = true;
    }
  }
}

double calculate_mean(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0f;
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0f);
    return sum / data.size();
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

std::string vec_to_pylist(const std::vector<uint32_t>& v){
    std::ostringstream ss;
    ss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i + 1 < v.size()) ss << ", ";
    }
    ss << "]";
    return ss.str();
}

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

void write_solution_status(std::string prefix, std::string solution_status){
  std::ofstream out(prefix + "_solution_status.tsv", std::ios::app);
  out << "solution_status\n";
  out << solution_status << "\n";
  out.close();
}

void write_variant_assignments(std::string prefix, std::vector<variant> variants){
  std::ofstream out(prefix + "_variant_assignments.tsv", std::ios::app);
  out << "POS\tNUC\tGAPPED_FREQ\tPOSTERIOR_PROB\tCOMPONENT_ASSIGMENT\tCONSENSUS_ASSIGNMENT\n";
  for(auto var : variants){
    if(var.assigned_component == -1) continue;
    out << std::to_string(var.position) << "\t";
    out << var.nuc << "\t";
    out << std::to_string(var.gapped_freq) << "\t";
    out << vec_to_pylist(var.marginal_posterior_probabilities) << "\t";
    out << std::to_string(var.assigned_component) << "\t";
    out << vec_to_pylist(var.consensus_numbers) << "\n";
  }
  out.close();
  
}

void write_solutions(std::string prefix, std::vector<std::vector<double>> solutions){
  std::ofstream out(prefix + "_solutions.tsv", std::ios::app);
  out << "solutions\n";
  for(auto solution : solutions){
    out << vec_to_pylist(solution) << "\n";
  }
  out.close(); 
}

void write_means(std::string prefix, std::vector<double> means){
  std::ofstream out(prefix + "_means.tsv", std::ios::app);
  out << "n\tmeans\n";
  out << std::to_string(means.size()) << "\t";
  out << vec_to_pylist(means) << "\n";
  out.close(); 
}

void write_error_limits(std::string prefix, double lower_bound, double upper_bound){
  std::ofstream out(prefix + "_error_limits.tsv", std::ios::app);
  out << "lower_error_bound\tupper_error_bound\n";
  out << std::to_string(lower_bound) << "\t";
  out << std::to_string(upper_bound) << "\n";
  out.close(); 
}

int gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, \
                              std::string ref, double default_threshold){
  min_qual = 9;
  if(ref.empty()){
    std::cerr << "Please provide a reference sequence." << std::endl;
    return(1);
  }
  double error_rate;
  std::vector<double> solution;
  std::vector<double> means;

  uint32_t round_val = 4;
  std::vector<variant> base_variants;
  parse_internal_variants(prefix, base_variants, min_depth, round_val, min_qual);
  set_deletion_flags(base_variants, 0.001);

  cluster_error(base_variants, min_qual, min_depth, error_rate);
  double lower_bound = 1-error_rate+0.0001;
  double upper_bound = error_rate-0.0001;
  write_error_limits(output_prefix, lower_bound, upper_bound);
  std::cerr << "upper bound " << upper_bound << " lower bound " << lower_bound << std::endl;
  set_freq_range_flags(base_variants, lower_bound, upper_bound, true);
  set_deletion_flags(base_variants, lower_bound);
  set_insertion_flags(base_variants);
  std::cerr << "filename: " << output_prefix << std::endl;
  std::cerr << "lower bound " << lower_bound <<  " upper bound " << upper_bound << std::endl;

  set_freq_range_flags(base_variants, lower_bound, upper_bound, true);

  std::vector<double> frequencies;
  std::vector<double> x_logit;
  std::vector<uint32_t> sites;

  for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].depth_flag) continue;
    if(base_variants[i].outside_freq_range) continue;
    if(base_variants[i].qual_flag) continue;
    std::cerr << base_variants[i].gapped_freq << " " << base_variants[i].position << " " << base_variants[i].qual << " " << base_variants[i].gapped_depth << std::endl;
    //if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked && base_variants[i].include_clustering){
      frequencies.push_back(base_variants[i].gapped_freq);
      sites.push_back(base_variants[i].position);
  }

  //handle the case of no variants less than the universal cluster
  if(frequencies.size() < 1){
    call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
    write_solution_status(output_prefix, "no solution:insufficient data to fit model");
    return(0);
  }

  int n_min = 0;
  std::unordered_map<int, int> counts;
  for (int v : sites) {
      ++counts[v];
  }

  std::unordered_map<int, int>::const_iterator it =
      std::max_element(
          counts.begin(), counts.end(),
          [](const std::pair<const int, int>& a,
            const std::pair<const int, int>& b) {
              return a.second < b.second;
          }
      );

  if (it != counts.end()) {
      n_min = it->second;
  }
  uint32_t N;
  if (n_min == 3){
    N = 6;
  } else if (n_min == 4) {
    N = 9;
  } else if (n_min == 2){
    N = 4;
  } else {
    std::cerr << "Insufficient data to fit GMM\n";
    write_solution_status(output_prefix, "no solution:insufficient data to fit model");
    return(1);
  }
  std::cerr << frequencies.size() << " variants to cluster." << std::endl;
  std::cerr << "Fitting " << N << " components to data." << std::endl;
  const double eps = 1e-6;
  std::vector<double> logL_history;
  gmm_1d::logit_transform(frequencies, x_logit, eps);
  //fit the model
  gmm_1d model(N);
  int return_code = model.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6,
      false,
      false
  );
  model.get_distinct_components_count(sites);
  uint32_t final_N = model.merged_means.size();

  //TESTLINES PRINTING
  const auto& m = model.get_means();
  gmm_1d::sigmoid_transform(m, means, eps);

  for(auto m : means){
    std::cerr << "mean " << m << std::endl;
  }

  std::cerr << "Final number of components after merging: " << final_N << std::endl;
  //fit the model a second time
  logL_history.clear();
  gmm_1d model2(final_N);
  int return_code2 = model2.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6,
      false,
      false
  );

  if(return_code2 == 1){
    final_N++;
    gmm_1d model2(final_N);
    int return_code2 = model2.fit(
        x_logit,
        sites,
        logL_history,
        20,
        1e-6,
        false,
        false
    );
  }

  means.clear();
  const auto& m2 = model2.get_means();
  gmm_1d::sigmoid_transform(m2, means, eps);

  for(auto m : means){
    std::cerr << "mean " << m << std::endl;
  }

  write_means(output_prefix, means);

  //attempt to solve the subset sum issue
  std::vector<std::vector<double>> solution_sets;
  bool solved = subset_sum(means, solution_sets); 
  if(!solved){
    std::cerr << "Could not solve unit sum problem for model means." << std::endl;
    call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
    write_solution_status(output_prefix, "no solution:unit-sum failed");
    return(1);
  }
  if(solution_sets.size() > 1){
    std::cerr << "Multiple solutions found for unit sum problem." << std::endl;
    call_majority_consensus(base_variants, output_prefix, default_threshold, min_depth);
    write_solution_status(output_prefix, "no solution:multiple possible solutions");
    write_solutions(output_prefix, solution_sets);
    return(1);
  }
  solution = solution_sets[0];
  
  //predict on all signals
  frequencies.clear();
  sites.clear();
  x_logit.clear();
  for(auto var : base_variants){
    if(var.depth_flag) continue;
    if(var.outside_freq_range) continue;
    if(var.qual_flag) continue;
    frequencies.push_back(var.gapped_freq);
    sites.push_back(var.position);
  }
  gmm_1d::logit_transform(frequencies, x_logit, eps);
  std::vector<uint32_t> assigned_components;
  std::vector<std::vector<double>> marginal_posterior_probabilities;
  model2.predict(
      x_logit,
      sites,
      assigned_components,
      marginal_posterior_probabilities
  );
  //add the assigments into the variants object
  uint32_t j = 0;
  for(uint32_t i=0; i < base_variants.size(); i++){
    variant &var = base_variants[i];
    if(var.depth_flag) continue;
    if(var.outside_freq_range) continue;
    if(var.qual_flag) continue;
    var.marginal_posterior_probabilities = marginal_posterior_probabilities[j];
    var.assigned_component = assigned_components[j];
    j++;
  }

  //use the solution set to assign consensus numbers
  std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> mapping_combinations;
  solve_additive_peaks(solution, means, mapping_combinations);
  assign_consensus_numbers(base_variants, mapping_combinations); 
  //check to make sure the cluster assignment is clear
  compare_component_assigments(base_variants);
  //place the universal cluster values in every consensus
  assign_universal_components(base_variants, upper_bound, means.size());
  cluster_consensus(base_variants, 
                    output_prefix, 
                    default_threshold,
                    min_depth, min_qual,
                    solution, means, 
                    ref, error_rate);
  write_solution_status(output_prefix, "solution");
  write_solutions(output_prefix, solution_sets);
  write_variant_assignments(output_prefix, base_variants);

  return 0;
}
