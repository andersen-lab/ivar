#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"

double compute_noise_threshold(const std::vector<double>& noise,
                               double m_min = 1,
                               double m_max = 10.0,
                               double s0 = 0.01,
                               double gamma = 0.5) {
    if (noise.empty()) return 0.0;

    // Compute mean
    double mean = std::accumulate(noise.begin(), noise.end(), 0.0) / noise.size();

    // Compute stddev
    double sq_sum = 0.0;
    for (double v : noise) {
        sq_sum += (v - mean) * (v - mean);
    }
    double sigma = std::sqrt(sq_sum / (noise.size() - 1));

    // Reference scale: if not given, use sigma itself
    if (s0 <= 0.0) s0 = sigma;

    // Adaptive multiplier: increases with sigma, saturates at m_max
    double m = m_min + (m_max - m_min) * std::pow(sigma / (sigma + s0), gamma);

    // Final cutoff depends *only* on sigma (noise level), not mean
    return m * sigma;
}


std::vector<double> z_score(std::vector<double> data) {
    double mean = calculate_mean(data);
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / data.size() - mean * mean);

    std::vector<double> z_scores;
    for (double x : data)
        z_scores.push_back((x - mean) / stddev);
    return z_scores;
}

std::vector<uint32_t>determine_outlier_points(std::vector<double> cluster){
    std::vector<uint32_t> removal_points;
    //calculate cluster specific percentiles
    std::vector<double> z_scores = z_score(cluster);
    for(uint32_t i=0; i < z_scores.size(); i++){
      double abs = std::abs(z_scores[i]);
      if(cluster[i] < 0.90) std::cerr << abs << " " << cluster[i] << std::endl;
      if(abs >= 5){
        std::cerr << "remove " << abs << " " << cluster[i] << std::endl;
        removal_points.push_back(i);
      }
    }
    return(removal_points);
}

void cluster_error(std::vector<variant> base_variants, uint8_t quality_threshold, uint32_t depth_cutoff, double &error_rate, double &error_std){
  double lower_bound = 0.50;
  double upper_bound = 0.99;
  set_freq_range_flags(base_variants, lower_bound, upper_bound, false);
  std::vector<variant> variants_original;
  uint32_t useful_count_original = 0;
  uint32_t max_pos = 0;
  std::vector<double> frequencies;

   for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].position > max_pos) max_pos = base_variants[i].position;
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
  uint32_t count_original=0;
  for(uint32_t i = 0; i < variants_original.size(); i++){
    double tmp = static_cast<double>(variants_original[i].gapped_freq);
    data_original.col(count_original) = tmp;
    count_original += 1;
  }
  //start with a small n value and if we don't find two major clusters we increase the number of clusters
  uint32_t n = 1;
  kmeans_model model;
  uint32_t chosen_peak = 0;

  while(n <= 2){
    model = train_model(n, data_original, true);
    std::vector<double> means = model.means;
    //index of largest mean
    uint32_t index = std::distance(means.begin(), std::max_element(means.begin(), means.end()));
    chosen_peak = index;
    bool stop=false;
    for(uint32_t i=0; i < model.clusters.size(); i++){
      double mean = calculate_mean(model.clusters[i]);
      double mad = calculate_mad(model.clusters[i], mean);
      double std = calculate_standard_deviation(model.clusters[i]);
      if(mad < 0.01 && i == index) {
        stop = true;
      }
      std::cerr << "n " << n << " mean " << mean << " mad " << mad << " " << model.clusters[i].size() << " std " << std << std::endl;
    }
    if(stop) break;
    else n++;
  }

  error_std = compute_noise_threshold(model.clusters[chosen_peak]);
  std::cerr << "error_std " << error_std  << std::endl;
  //exit(0);
  //for each cluster this describes the points which are outliers
  std::vector<uint32_t> outliers = determine_outlier_points(model.clusters[chosen_peak]);
  //std::vector<uint32_t> outliers;
  std::vector<double> universal_cluster = model.clusters[chosen_peak];
  std::vector<double> cleaned_cluster;
  for(uint32_t i=0; i < universal_cluster.size(); i++){
    auto it = std::find(outliers.begin(), outliers.end(), i);
    if(it == outliers.end()){
      cleaned_cluster.push_back(universal_cluster[i]);
    }
  }

  //get the upper edge of the noise cluster
  auto min_it = std::min_element(cleaned_cluster.begin(), cleaned_cluster.end());
  error_rate = *min_it;
}
