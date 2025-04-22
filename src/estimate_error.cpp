#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"

std::vector<double> z_score(std::vector<double> data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / data.size() - mean * mean);

    std::vector<double> z_scores;
    for (double x : data)
        z_scores.push_back((x - mean) / stddev);
    return z_scores;
}

std::vector<std::vector<uint32_t>> determine_outlier_points(std::vector<double> data, std::vector<double> means){
    std::vector<std::vector<double>> clusters(means.size());
    //assign points to clusters
    for(auto tmp : data){
      auto it = std::min_element(means.begin(), means.end(), 
          [tmp](double a, double b) {
              return std::fabs(a - tmp) < std::fabs(b - tmp);
          });
      uint32_t index = std::distance(means.begin(), it);
      clusters[index].push_back(tmp);
    }

    std::vector<std::vector<uint32_t>> removal_points(means.size());
    //calculate cluster specific percentiles
    for(uint32_t j=0; j < clusters.size(); j++){
      std::vector<double> z_scores = z_score(clusters[j]);
      for(uint32_t i=0; i < z_scores.size(); i++){
        double abs = std::abs(z_scores[i]);
        if(abs >= 10){
          removal_points[j].push_back(i);
        }
      }
    }
    return(removal_points);
}

double cluster_error(std::string filename, uint8_t quality_threshold, uint32_t depth_cutoff){
  /*
    Here we use clustering to determine the value of the noise.
  */

  float lower_bound = 0.50;
  float upper_bound = 0.99;
  uint32_t round_val = 4;
  
  std::vector<variant> base_variants;
  parse_internal_variants(filename, base_variants, depth_cutoff, lower_bound, upper_bound, round_val, quality_threshold);
  if(base_variants[0].version_1_var){
    calculate_reference_frequency(base_variants, filename, depth_cutoff, lower_bound, upper_bound);
  }

  std::vector<variant> variants_original;
  uint32_t useful_count_original = 0;
  uint32_t max_pos = 0;
  std::vector<double> frequencies;
   for(uint32_t i=0; i < base_variants.size(); i++){
    if(base_variants[i].position > max_pos) max_pos = base_variants[i].position;
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && !base_variants[i].amplicon_masked && !base_variants[i].primer_masked){      
      useful_count_original++;
      variants_original.push_back(base_variants[i]);
      frequencies.push_back(base_variants[i].freq);
    }
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

  while(n <= 5){
    model = train_model(n, data_original);
    std::vector<double> means = model.means;
    bool stop=true;
    for(uint32_t i=0; i < model.clusters.size(); i++){
      float percent = (float)model.clusters[i].size() / (float)data_original.size();
      if(percent > 0.85) {
        n++;
        stop = false;
      }
      //std::cerr << "n " << n << " percent " << percent << std::endl;
    }
    if(stop) break;
  }
  //for each cluster this describes the points which are outliers
  std::vector<std::vector<uint32_t>> removal_points = determine_outlier_points(frequencies, model.means);

  uint32_t j = 0;
  uint32_t largest=0;
  for(uint32_t i=0; i < model.clusters.size(); i++){
    //std::cerr << "percent " << (float)model.clusters[i].size() / (float)data_original.size() << std::endl;
    //double mean = std::accumulate(model.clusters[i].begin(), model.clusters[i].end(), 0.0) / model.clusters[i].size();
    //std::cerr << "mean " << mean << std::endl;
    //double stdev = calculate_standard_deviation(model.clusters[i]);
    //std::cerr << stdev << std::endl;
    if(model.clusters[i].size() > largest){
      j = i;
      largest = model.clusters[i].size();
    }
  }
  std::vector<double> universal_cluster = model.clusters[j];
  std::vector<uint32_t> outliers = removal_points[j];
  std::vector<double> cleaned_cluster;
  for(uint32_t i=0; i < universal_cluster.size(); i++){
    auto it = std::find(outliers.begin(), outliers.end(), i);
    if(it == outliers.end()){
      cleaned_cluster.push_back(universal_cluster[i]);
    }
  }

  
  //get the upper edge of the noise cluster
  auto min_it = std::min_element(cleaned_cluster.begin(), cleaned_cluster.end());
  return(*min_it);
}
