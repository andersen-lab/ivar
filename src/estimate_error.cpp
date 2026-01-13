#include "estimate_error.h"
#include "gmm.h"
#include "saga.h"
#include "gmm_1d.h"
#include <numeric>
std::vector<double> z_score(std::vector<double> data) {
    double mean = calculate_mean(data);
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / data.size() - mean * mean);
    std::vector<double> z_scores;
    for (double x : data)
        z_scores.push_back((x - mean) / stddev);
    return z_scores;
}

void determine_outlier_points(std::vector<double> cluster, double threshold, std::vector<uint32_t> &removal_points){
  std::vector<double> z_scores = z_score(cluster);
  for(uint32_t i=0; i < z_scores.size(); i++){
    double abs = std::abs(z_scores[i]);
    if(abs >= threshold){
      std::cerr << "remove " << abs << " " << cluster[i] << std::endl;
      removal_points.push_back(i);
    }
  }
}

void cluster_error(std::vector<variant> base_variants, uint8_t quality_threshold, uint32_t depth_cutoff, double &error_rate){
  double lower_bound = 0.50;
  double upper_bound = 0.99;
  const double eps = 1e-6;
  set_freq_range_flags(base_variants, lower_bound, upper_bound, false);

  std::vector<double> frequencies;
  std::vector<uint32_t> sites;

  for(uint32_t i=0; i < base_variants.size(); i++){
    if(!base_variants[i].amplicon_flux && !base_variants[i].depth_flag && \
      !base_variants[i].outside_freq_range && !base_variants[i].qual_flag && \
      !base_variants[i].amplicon_masked){
      frequencies.push_back(base_variants[i].gapped_freq);
      sites.push_back(base_variants[i].position);
    }
  }

  if(frequencies.empty()){
    error_rate = 0.01;
    return;
  }

  std::vector<double> x_logit;
  gmm_1d::logit_transform(frequencies, x_logit, eps);
  int N = 5;

  // Fit GMM
  std::vector<double> logL_history;
  gmm_1d model(N);
  model.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6,
      true,
      false
  );

  model.get_distinct_components_count(sites);

  std::cerr << "Merged components: " << model.merged_means.size() << "\n";
  
  uint32_t final_N = model.merged_means.size();
  gmm_1d model2(final_N);
  model2.fit(
      x_logit,
      sites,
      logL_history,
      20,
      1e-6,
      true
  );
  const auto& m2 = model2.get_means();
  std::vector<double> m_sigmoid2;

  gmm_1d::sigmoid_transform(m2, m_sigmoid2, eps);
  double largest_mean = 0.0;
  uint32_t largest_idx = 0;
  for(uint32_t m=0; m < m_sigmoid2.size(); ++m){
    std::cerr << "mean " << m_sigmoid2[m] << "\n";
    if(m_sigmoid2[m] > largest_mean){
      largest_mean = m_sigmoid2[m];
      largest_idx = m;
    }
  }
  std::vector<uint32_t> assigned_components;
  std::vector<std::vector<double>> marginal_posterior_probabilities;
  model2.predict(
      x_logit,
      sites,
      assigned_components,
      marginal_posterior_probabilities
  );

  std::vector<double> largest_cluster_snps;
  for(uint32_t i =0; i < frequencies.size(); i++){
    if(assigned_components[i] == largest_idx){
      largest_cluster_snps.push_back(frequencies[i]);
    }
  }

  double smallest_cluster_snp = 1.0;
  for(uint32_t k=0; k < largest_cluster_snps.size(); k++){
    if(largest_cluster_snps[k] < smallest_cluster_snp){
      smallest_cluster_snp = largest_cluster_snps[k];
    }
  }
  error_rate = smallest_cluster_snp;
}
