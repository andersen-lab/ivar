#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"


// Map input x [0.001,0.10] to output y [1,30], decreasing nonlinearly
double adaptive_value(double x,
                      double y_min = 1.0,
                      double y_max = 15.0,
                      double x0 = 0.015,
                      double gamma = 1.5) {

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

std::vector<uint32_t>determine_outlier_points(std::vector<double> cluster, double threshold){
    std::vector<uint32_t> removal_points;
    //calculate cluster specific percentiles
    std::vector<double> z_scores = z_score(cluster);
    for(uint32_t i=0; i < z_scores.size(); i++){
      double abs = std::abs(z_scores[i]);
      if(abs >= threshold){
        //std::cerr << "remove " << abs << " " << cluster[i] << std::endl;
        removal_points.push_back(i);
      }
    }
    return(removal_points);
}

bool is_within_std(const std::vector<double>& values, double mean, double std_dev) {
  for (double val : values) {
    if (std::abs(val - mean) < std_dev) {
      return true;  //found at least one value within the standard deviation
    }
  }
  return false; //no value within the standard deviation
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

  uint32_t n = 5;
  kmeans_model model;
  uint32_t chosen_peak = 0;
  bool solved = true;
  while(n >= 2){
    bool close_clusters = false;
    solved = true;
    model = train_model(n, data_original, true);
    std::vector<double> means = model.means;
    //index of largest cluster
    chosen_peak = std::distance(means.begin(), std::max_element(means.begin(), means.end()));

    std::vector<uint32_t> indices(means.size());
    for (uint32_t i = 0; i < means.size(); ++i) {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](uint32_t a, uint32_t b) {
      return means[a] > means[b];
    });
    for(auto m : means){
      std::cerr << "m " << m << std::endl;
    }
    //compute centroids and within-cluster standard deviation
    for (size_t i = 0; i < 2; ++i) {
      uint32_t index = indices[i];
      double centroid = means[index];
      //remove this value from the means
      std::vector<double> tmp_means = means;
      if (index < tmp_means.size()) {
        tmp_means.erase(tmp_means.begin() + index);
      }

      //make sure there isn't another cluster within a couple standard deviations
      double within_std = calculate_standard_deviation(model.clusters[index]);
      close_clusters = is_within_std(tmp_means,  centroid, within_std*3);
      if(close_clusters){
        solved = false;
        std::cerr << "here" << std::endl;
        break;
      }
    }
    if(solved){
      break;
    }
    std::cerr << "here" << std::endl;
    n--;
  }
  std::cerr << "chosen n " << n << " chosen peak " << chosen_peak << std::endl;
  //exit(0);
  std::vector<double> cleaned_cluster;

  //for each cluster this describes the points which are outliers
  if(n > 1){
    std::vector<uint32_t> outliers = determine_outlier_points(model.clusters[chosen_peak], 2);
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
    //TODO handle the case of n=1
    std::vector<uint32_t> outliers = determine_outlier_points(frequencies, 2);
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
  std::cerr << *min_it << std::endl;
  error_rate = *min_it;
}
