#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"

std::vector<double> z_score(std::vector<double> data) {
    double mean = calculate_mean(data);
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / data.size() - mean * mean);
    std::vector<double> z_scores;
    for (double x : data)
        z_scores.push_back((x - mean) / stddev);
    return z_scores;
}

std::vector<uint32_t>determine_outlier_points(std::vector<double> cluster, double threshold){
    std::vector<uint32_t> removal_points;
    std::vector<double> z_scores = z_score(cluster);
    for(uint32_t i=0; i < z_scores.size(); i++){
      double abs = std::abs(z_scores[i]);
      if(abs >= threshold){
        std::cerr << "remove " << abs << " " << cluster[i] << std::endl;
        removal_points.push_back(i);
      }
    }
    return(removal_points);
}
void cluster_error(std::vector<variant> base_variants, uint8_t quality_threshold, uint32_t depth_cutoff, double &error_rate){
  double lower_bound = 0.50;
  double upper_bound = 0.99;
  set_freq_range_flags(base_variants, lower_bound, upper_bound, false);

  std::vector<variant> variants_original;
  uint32_t useful_count_original = 0;
  std::vector<double> frequencies;
  std::vector<uint32_t> subsample_position;

  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && \
      !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && \
      !base_variants[i].amplicon_masked){
      useful_count_original++;
      variants_original.push_back(base_variants[i]);
      frequencies.push_back(base_variants[i].gapped_freq);
      subsample_position.push_back(base_variants[i].position);
    }
  }

  if(variants_original.empty()){
    error_rate = 1;
    return;
  }
  std::cerr << "error estimate useful var " << useful_count_original << std::endl;
  arma::mat data_original(1, useful_count_original, arma::fill::zeros);
  uint32_t count_original=0;
  for(uint32_t i = 0; i < variants_original.size(); i++){
    double tmp = static_cast<double>(variants_original[i].gapped_freq);
    data_original.col(count_original) = tmp;
    count_original += 1;
  }

  uint32_t chosen_peak = 0;
  double var_floor = 0.0001;
  bool clustering_failed = false;
  std::vector<double> all_bics;
  uint32_t bootstrap_reps = 1000;
  uint32_t i = 1;
  std::vector<variant> subsampled_variants;
  std::unordered_map<uint32_t, uint32_t> model_counter; 
  while(i <= bootstrap_reps){
    uint32_t n = 1;
    subsampled_variants.clear();
    clustering_failed = false;
    all_bics.clear();
    arma::mat subsample = subsample_with_replacement(data_original, data_original.size(), subsample_position, subsampled_variants, variants_original);
    while(n <= 5){
      reset_variants_info(subsampled_variants);
      gaussian_mixture_model model = retrain_model(n, subsample, subsampled_variants, 2, var_floor, clustering_failed);
      for(auto cluster : model.clusters){
        if(cluster.size() == 0){
          clustering_failed = true;
        }
      }
      if(!clustering_failed){
        all_bics.push_back(model.bic);
      } else {
        all_bics.push_back(std::numeric_limits<double>::max());
      }
      n++;
    }
    double lowest = std::numeric_limits<double>::max(); 
    uint32_t winner;
    for(uint32_t i=0; i < all_bics.size(); i++){
      if(all_bics[i] < lowest){
        lowest = all_bics[i];
        winner = i+1;
      }
    }
    model_counter[winner] += 1;
    i++;
  }

  uint32_t optimal_n;
  uint32_t highest = 0;
  for (const auto& [key, value] : model_counter) {
    if(value >= highest){
      optimal_n = key;
      highest = value;
    }
    std::cerr << key << " : " << value << '\n';
  }

  gaussian_mixture_model model = retrain_model(optimal_n, data_original, variants_original, 2, var_floor, clustering_failed);
  std::vector<double> means = model.means;
  for(auto m : means){
    std::cerr << "mean " << m << std::endl;
  }
  chosen_peak = std::distance(means.begin(), std::max_element(means.begin(), means.end()));
  std::vector<double> cleaned_cluster;
  std::vector<uint32_t> outliers;

  //for each cluster this describes the points which are outliers
  if(optimal_n > 1){
    outliers = determine_outlier_points(model.clusters[chosen_peak], 2.5);
    std::vector<double> universal_cluster = model.clusters[chosen_peak];
    for(uint32_t i=0; i < universal_cluster.size(); i++){
      auto it = std::find(outliers.begin(), outliers.end(), i);
      if(it == outliers.end()){
        cleaned_cluster.push_back(universal_cluster[i]);
      }
    }
  } else {
    //TODO handle the case of n=1
    outliers = determine_outlier_points(frequencies, 2.5);
    for(uint32_t i=0; i < frequencies.size(); i++){
      auto it = std::find(outliers.begin(), outliers.end(), i);
      if(it == outliers.end()){
        cleaned_cluster.push_back(frequencies[i]);
      }
    }
  }

  auto min_it = std::min_element(cleaned_cluster.begin(), cleaned_cluster.end());
  error_rate = *min_it;
}
