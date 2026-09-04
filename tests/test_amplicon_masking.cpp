#include <iostream>
#include <random>
#include <vector>

#include "../src/gmm.h"
#include "../src/gmm_1d.h"

void amplicon_specific_cluster_assignment(std::vector<variant> &variants, const gmm_1d &model, const std::vector<int> &component_indices);

// Pinned explicitly rather than relying on DEFAULT_AMPLICON_STDEV: that default
// is 2.0 (masking off), so these cases only exercise anything with a real
// threshold passed in.
static const double TEST_STDEV = 0.05;

variant make_variant(uint32_t position, const std::vector<std::string> &ids,
                     const std::vector<double> &freqs, const std::vector<uint32_t> &depths) {
  variant v;
  v.position = position;
  v.nuc = "A";
  v.amplicon_ids = ids;
  v.freq_numbers = freqs;
  v.amplicon_depths = depths;
  return v;
}

bool expect(bool condition, const std::string &message) {
  if(!condition) std::cerr << "FAIL: " << message << std::endl;
  return condition;
}

// The stdev flag is a prefilter on the per-amplicon frequencies
bool test_stdev_flagging() {
  std::vector<variant> variants = {
    make_variant(10, {"1-500", "300-800"}, {0.50, 0.52}, {100, 100}),
    make_variant(20, {"1-500", "300-800"}, {0.50, 0.10}, {100, 100}),
    make_variant(30, {"1-500"}, {}, {}),
  };
  flag_amplicon_variation(variants, TEST_STDEV);

  bool ok = true;
  ok &= expect(!variants[0].position_masked, "tight frequencies should not be masked");
  ok &= expect(variants[1].position_masked, "split frequencies should be masked");
  ok &= expect(!variants[2].position_masked, "a single amplicon cannot show flux");
  return ok;
}

// Identical frequencies must be able to land either side of the threshold purely on
// the strength of their depths
bool test_depth_weighting() {
  std::vector<variant> even = { make_variant(10, {"1-500", "300-800"}, {0.50, 0.62}, {100, 100}) };
  std::vector<variant> lopsided = { make_variant(10, {"1-500", "300-800"}, {0.50, 0.62}, {1000, 10}) };
  flag_amplicon_variation(even, TEST_STDEV);
  flag_amplicon_variation(lopsided, TEST_STDEV);

  bool ok = true;
  ok &= expect(even[0].position_masked, "evenly weighted split should be masked");
  ok &= expect(!lopsided[0].position_masked, "a shallow outlier amplicon should be down-weighted");
  return ok;
}

// Masking must not depend on the order variants appear in, which is what the
// write-time implementation could not guarantee
bool test_masking_is_order_independent() {
  variant low = make_variant(100, {"1-500"}, {}, {});
  variant high = make_variant(400, {"1-500", "300-800"}, {0.50, 0.10}, {100, 100});

  std::vector<variant> forward = {low, high};
  std::vector<variant> reverse = {high, low};
  flag_amplicon_variation(forward, TEST_STDEV);
  flag_amplicon_variation(reverse, TEST_STDEV);

  bool ok = true;
  ok &= expect(forward[0].amplicon_masked, "flux at a later position must mask an earlier one");
  ok &= expect(reverse[1].amplicon_masked, "flux at an earlier position must mask a later one");
  ok &= expect(forward[0].amplicon_masked == reverse[1].amplicon_masked,
               "masking must not depend on variant order");
  return ok;
}

// amplicon_masked reports flux found elsewhere, so a variant's own flux must not set it
bool test_self_flux_does_not_mask() {
  std::vector<variant> alone = { make_variant(400, {"1-500", "300-800"}, {0.50, 0.10}, {100, 100}) };
  flag_amplicon_variation(alone, TEST_STDEV);

  std::vector<variant> pair = {
    make_variant(400, {"1-500", "300-800"}, {0.50, 0.10}, {100, 100}),
    make_variant(450, {"1-500"}, {0.50, 0.10}, {100, 100}),
  };
  flag_amplicon_variation(pair, TEST_STDEV);

  bool ok = true;
  ok &= expect(alone[0].position_masked, "the lone variant should still be position masked");
  ok &= expect(!alone[0].amplicon_masked, "a variant must not be amplicon masked by its own flux");
  ok &= expect(pair[0].amplicon_masked, "a second flagged variant on a shared amplicon should mask");
  return ok;
}

gmm_1d fit_bimodal_model(std::vector<int> &component_indices) {
  std::mt19937 rng(112358);
  std::normal_distribution<double> low(0.20, 0.005);
  std::normal_distribution<double> high(0.80, 0.005);

  std::vector<double> x;
  for(int i = 0; i < 500; i++){
    x.push_back(low(rng));
    x.push_back(high(rng));
  }

  gmm_1d model(12, 42);
  model.set_use_half_normal_for_noise(true, 0.97);
  model.set_mean_precision_prior(0.5);
  model.fit(x);
  model.set_min_cluster_fraction(0.1);
  component_indices = model.get_effective_components(model.predict(x));
  return model;
}

// The GMM check clears the stdev flag when the per-amplicon frequencies all describe
// the same population, and leaves it when they genuinely span clusters
bool test_cluster_agreement_rewrites_flag() {
  std::vector<int> component_indices;
  gmm_1d model = fit_bimodal_model(component_indices);
  if(component_indices.size() < 2){
    std::cerr << "FAIL: fixture model did not resolve two components" << std::endl;
    return false;
  }

  std::vector<variant> variants = {
    make_variant(10, {"1-500", "300-800"}, {0.19, 0.21}, {100, 100}),
    make_variant(20, {"1-500", "300-800"}, {0.20, 0.80}, {100, 100}),
  };
  for(auto &v : variants) v.position_masked = true;

  amplicon_specific_cluster_assignment(variants, model, component_indices);
  rewrite_position_masking(variants);

  bool ok = true;
  ok &= expect(!variants[0].position_masked, "one cluster across amplicons is not flux");
  ok &= expect(variants[1].position_masked, "two clusters across amplicons is flux");
  return ok;
}

// Assignments are argmaxed over the components the model kept. Ranging over all
// components lets one the model discarded split two frequencies of a single cluster:
// in this fixture 0.35 argmaxes onto a discarded component near 0.21, while both 0.20
// and 0.35 belong to the same effective component
bool test_assignments_use_effective_components() {
  std::vector<int> component_indices;
  gmm_1d model = fit_bimodal_model(component_indices);

  std::vector<variant> variants = { make_variant(10, {"1-500", "300-800"}, {0.20, 0.35}, {100, 100}) };
  flag_amplicon_variation(variants, TEST_STDEV);
  if(!expect(variants[0].position_masked, "fixture should be flagged by stdev first")) return false;

  amplicon_specific_cluster_assignment(variants, model, component_indices);
  rewrite_position_masking(variants);

  bool ok = true;
  for(uint32_t assignment : variants[0].freq_assignments){
    bool effective = false;
    for(int ci : component_indices){
      if((int)assignment == ci) effective = true;
    }
    ok &= expect(effective, "assignment fell outside the effective components");
  }
  ok &= expect(!variants[0].position_masked,
               "a discarded component must not split one cluster into two");
  return ok;
}

// Assignments are rebuilt rather than appended to
bool test_assignments_do_not_accumulate() {
  std::vector<int> component_indices;
  gmm_1d model = fit_bimodal_model(component_indices);

  std::vector<variant> variants = { make_variant(10, {"1-500", "300-800"}, {0.20, 0.80}, {100, 100}) };
  variants[0].position_masked = true;

  amplicon_specific_cluster_assignment(variants, model, component_indices);
  size_t once = variants[0].freq_assignments.size();
  amplicon_specific_cluster_assignment(variants, model, component_indices);

  return expect(variants[0].freq_assignments.size() == once,
                "repeated assignment should not accumulate entries");
}

int main() {
  bool ok = true;
  ok &= test_stdev_flagging();
  ok &= test_depth_weighting();
  ok &= test_masking_is_order_independent();
  ok &= test_self_flux_does_not_mask();
  ok &= test_cluster_agreement_rewrites_flag();
  ok &= test_assignments_use_effective_components();
  ok &= test_assignments_do_not_accumulate();
  return ok ? 0 : -1;
}
