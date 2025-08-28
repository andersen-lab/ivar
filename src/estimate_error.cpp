#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"


// Map input x [0.001,0.10] to output y [1,30], decreasing nonlinearly
double adaptive_value(double x,
                      double y_min = 1.0,
                      double y_max = 15.0,
                      double x0 = 0.015,
                      double gamma = 2.0) {

    // Power-law / inverse mapping
    double y = y_min + (y_max - y_min) / (std::pow(x / x0, gamma) + 1.0);
    return y*x;
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
      if(abs >= 3){
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
  uint32_t n = 2;
  kmeans_model model;
  uint32_t chosen_peak = 0;


  double std = calculate_standard_deviation(frequencies);
  uint32_t real_n = 0;
  double err = 0;
  while(n <= 2){
    model = train_model(n, data_original, true);
    std::vector<double> means = model.means;

    double mean_diff = std::abs(means[0] - means[1]);
    if(mean_diff > std){
      real_n = 2;
    } else {
      real_n = 1;
    }
    //index of largest mean
    uint32_t index = std::distance(means.begin(), std::max_element(means.begin(), means.end()));
    chosen_peak = index;
    for(uint32_t i=0; i < means.size(); i++){
      double std_1 = calculate_standard_deviation(model.clusters[i]);
      std::cerr << means[i] << " " << std_1 << std::endl;
    }

    n++;
  }
  std::cerr << real_n << std::endl;
  std::vector<double> cleaned_cluster;
  //for each cluster this describes the points which are outliers
  if(real_n == 2){
    std::vector<uint32_t> outliers = determine_outlier_points(model.clusters[chosen_peak]);
    std::vector<double> universal_cluster = model.clusters[chosen_peak];
    for(uint32_t i=0; i < universal_cluster.size(); i++){
      auto it = std::find(outliers.begin(), outliers.end(), i);
      if(it == outliers.end()){
        cleaned_cluster.push_back(universal_cluster[i]);
      }
    }
    //for(auto f : cleaned_cluster) std::cerr << f << std::endl;
    double cstd = calculate_standard_deviation(cleaned_cluster);
    std::cerr << "c standard " << cstd<< std::endl;
    error_std = adaptive_value(cstd);
  } else {
    std::vector<uint32_t> outliers = determine_outlier_points(frequencies);
    for(uint32_t i=0; i < frequencies.size(); i++){
      auto it = std::find(outliers.begin(), outliers.end(), i);
      if(it == outliers.end()){
        cleaned_cluster.push_back(frequencies[i]);
      }
    }
    double cstd = calculate_standard_deviation(cleaned_cluster);
    std::cerr << "c standard " << cstd<< std::endl;
    error_std = adaptive_value(cstd);

  }
  auto min_it = std::min_element(cleaned_cluster.begin(), cleaned_cluster.end());
  error_rate = *min_it;
}
