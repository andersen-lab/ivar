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
    //calculate cluster specific percentiles
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

   for(uint32_t i=0; i < base_variants.size(); i++){

    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked){
      useful_count_original++;
      variants_original.push_back(base_variants[i]);
      frequencies.push_back(base_variants[i].gapped_freq);
    }
  }
  if(variants_original.empty()){
    error_rate = 1;
    return;
  }
  arma::mat data_original(1, useful_count_original, arma::fill::zeros);
  arma::mat data_second_dim(2, useful_count_original, arma::fill::zeros); 
  uint32_t count_original=0;
  for(uint32_t i = 0; i < variants_original.size(); i++){
    double tmp = static_cast<double>(variants_original[i].gapped_freq);
    data_original.col(count_original) = tmp;
    data_second_dim(0, i) = tmp;
    data_second_dim(1, i) = 0;
    count_original += 1;
  }

  uint32_t n = 1;
  gaussian_mixture_model model;
  uint32_t chosen_peak = 0;

  bool clustering_failed = false;
  uint32_t optimal_n=0;
  std::vector<double> bics;
  std::vector<double> ns;
  while(n <= 5){
    std::cerr << "error estimate " << n << std::endl;

    model = retrain_model(n, data_original, data_second_dim, variants_original, 2, 0.001, clustering_failed);
    if(clustering_failed) {
      n++;
      continue;
    }
    for(auto cluster : model.clusters){
      if(cluster.size() == 0){
        clustering_failed = true;
      }
    }
    if(clustering_failed) {
      n++;
      continue;
    }
    if(n == 1){
      double std = calculate_standard_deviation(model.clusters[0]);
      if(std < 0.03){
        optimal_n = 1;
        break;
      }
    }

    bics.push_back(model.bic);
    ns.push_back((double)n);
    std::cerr << "bic " << model.bic << std::endl;
    std::vector<double> means = model.means;
    for(uint32_t i=0; i < means.size(); i++){
      std::cerr << means[i] << " " << calculate_standard_deviation(model.clusters[i]) << std::endl;
    }
    //index of largest cluster
    chosen_peak = std::distance(means.begin(), std::max_element(means.begin(), means.end()));
    n++;
  }
  if(optimal_n == 0){
    std::vector<double> exclude_ns;
    optimal_n = elbow_method(bics, ns, exclude_ns);
  }
  std::cerr << "optimal n " << optimal_n << std::endl;

  model = retrain_model(optimal_n, data_original, data_second_dim, variants_original, 2, 0.001, clustering_failed);
  std::vector<double> means = model.means;
  chosen_peak = std::distance(means.begin(), std::max_element(means.begin(), means.end()));

  std::vector<double> cleaned_cluster;
  std::vector<uint32_t> outliers;

  //for each cluster this describes the points which are outliers
  if(n > 1){
    //outliers = determine_outlier_points(model.clusters[chosen_peak], 3);
    std::vector<double> universal_cluster = model.clusters[chosen_peak];
    for(uint32_t i=0; i < universal_cluster.size(); i++){
      auto it = std::find(outliers.begin(), outliers.end(), i);
      if(it == outliers.end()){
        cleaned_cluster.push_back(universal_cluster[i]);
      }
    }
  } else {
    //TODO handle the case of n=1
    outliers = determine_outlier_points(frequencies, 3);
    for(uint32_t i=0; i < frequencies.size(); i++){
      auto it = std::find(outliers.begin(), outliers.end(), i);
      if(it == outliers.end()){
        cleaned_cluster.push_back(frequencies[i]);
      }
    }
  }

  auto min_it = std::min_element(cleaned_cluster.begin(), cleaned_cluster.end());
  std::cerr << *min_it << std::endl;
  error_rate = *min_it;
}
