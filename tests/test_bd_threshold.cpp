#include <iostream>
#include <algorithm>
#include <fstream>

#include "../src/gaussian_sampler.h"
#include "../src/gmm_1d.h"

int main() {
  int ncases = 5;
//  double mean_diffs[] = {0.1, 0.2, 0.3, 0.4, 0.5};
//  int n_samples[] = {10, 50, 100, 500,  1000};
//  double mean1[] = {0.1, 0.2, 0.3, 0.4, 0.5};
//  double stdevs[] = {0.025, 0.05, 0.075, 0.1};
//  double stdevs_diff[] = {0, 0.05, 0.075, 0.1};

  double mean_diffs[] = {0.1, 0.2, 0.3, 0.4, 0.5};
  int n_samples[] = {100};
  double mean1[] = {0.1, 0.2, 0.3, 0.4, 0.5};
  double stdevs[] = { 0.05};
  double stdevs_diff[] = {0};

  std::ofstream f("/Users/karthik/tmp/test_bd_threshold/data/gaussian_samples_bd_threshold.tsv", std::ios::out);
  f << "mean1\tmean2\tstdev1\tstdev2\tnsamples\tsample\tvalue\n";

  std::ofstream fresults("/Users/karthik/tmp/test_bd_threshold/output/inferred_gmm.tsv", std::ios::out);
  fresults << "mean1\tmean2\tstdev1\tstdev2\tnsamples\tinferred_mean1\tinferred_mean2\tinferred_var1\tinferred_var2\tBD\tdistinct_components\n";

  for(int nsamples: n_samples) {
    // Set number of sites
    std::vector<int> sites;
    for(int j = 0;j < 2; j++)
      for(int i = 0;i < nsamples; i++)
        sites.push_back(i);

    for(double mean_diff : mean_diffs){
      for(double m1: mean1){
        for(double stdev1: stdevs) {
          for(double stdev_diff: stdevs_diff) {
            double stdev2 = stdev1 + stdev_diff;
            // Generate samples from two gaussian mixtures
            std::vector<double> samples1;
            std::vector<double> samples2;
            double m2 = m1 + mean_diff;
            gaussian_sampler sampler1(m1, stdev1, 42);
            gaussian_sampler sampler2(m2, stdev2, 42);
            sampler1.sample(samples1, nsamples);
            sampler2.sample(samples2, nsamples);

            // Clamp values to [0, 1]
            for(int n = 0; n < nsamples;n++){
              // Rudimentary fix to avoid values too close to 0s and 1s.
              if (samples1[n] < 0.0) samples1[n] = 1e-2;
              if (samples1[n] > 1.0) samples1[n] = 1 - 1e-2;

              if (samples2[n] < 0.0) samples2[n] = 1e-2;
              if (samples2[n] > 1.0) samples2[n] = 1-1e-2;
            }

            std::vector<double> samples;

            samples.insert(samples.end(), samples1.begin(), samples1.end());
            samples.insert(samples.end(), samples2.begin(), samples2.end());

            // Logit transform
//            std::vector<double> samples1_logit;
//            std::vector<double> samples2_logit;
//            gmm_1d::logit_transform(samples1, samples1_logit);
//            gmm_1d::logit_transform(samples2, samples2_logit);

//            std::vector<double> samples_logit;
//            samples_logit.insert(samples_logit.end(), samples1_logit.begin(), samples1_logit.end());
//            samples_logit.insert(samples_logit.end(), samples2_logit.begin(), samples2_logit.end());

            // Write to file
            for(int n = 0; n < nsamples;n++){
              f << m1 << "\t"
                << m2 << "\t"
                << stdev1 << "\t"
                << stdev2 << "\t"
                << nsamples << "\t"
                << "One" << "\t"
                << samples1[n] << "\n";
//                << samples1_logit[n] << "\n";
            }

            for(int n = 0; n < nsamples;n++){
              f << m1 << "\t"
                << m2 << "\t"
                << stdev1 << "\t"
                << stdev2 << "\t"
                << nsamples << "\t"
                << "Two" << "\t"
                << samples2[n] << "\n";
//                << samples2_logit[n] << "\n";
            }

            // Fit GMM
            std::vector<double> logL_history;
            gmm_1d model(2);
            model.set_use_half_normal_for_noise(false);
            model.fit(
                samples,
                sites,
                logL_history,
                {},
                20,
                1e-6
            );

            // Print results (logit space)
            const auto& w = model.get_weights();
            const auto& m = model.get_means();
            const auto& v = model.get_vars();

            std::vector<double> m_sigmoid;

            gmm_1d::sigmoid_transform(m, m_sigmoid);

            fresults << m1 << "\t"
                     << m2 << "\t"
                     << stdev1 << "\t"
                     << stdev2 << "\t"
                     << nsamples << "\t"
                     << m[0] << "\t"
                     << m[1] << "\t"
//                     << m_sigmoid[0] << "\t"
//                     << m_sigmoid[1] << "\t"
                     << v[0] << "\t"
                     << v[1] << "\t"
                     << gmm_1d::calculate_bhattacharyya_distance_1d(m[0], v[0], m[1], v[1]) << "\t"
                     << model.get_distinct_components_count(sites, 1.0) << "\n";
          }
        }
      }
    }

  }


  f.close();
  fresults.close();

  return 0;
}