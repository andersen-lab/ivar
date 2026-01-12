#include "site_state_aggregator_stats.h"
#include <cmath>

const double site_state_aggregator_stats::STDEV_THRESHOLD = 0.055;

void site_state_aggregator_stats::add_site(site_state &ss) {
  //TODO: Add explicit check if site_state.state == this->state
  site_amplicon_aggregator_stats *ssa_stats = find_or_create_site_amplicon_aggregator_stats(ss.amplicon);
  ssa_stats->add_site(ss.quality);
}

const site_amplicon_aggregator_stats* site_state_aggregator_stats::get_site_amplicon_aggregator_stats(ITNode* amplicon) const {
  for(auto &ait : site_state_amplicon_stats) {
    if (ait.get_amplicon() == amplicon) {
      return &ait;
    }
  }
  return nullptr;
}

site_amplicon_aggregator_stats* site_state_aggregator_stats::find_or_create_site_amplicon_aggregator_stats(ITNode *amplicon) {
  for(auto &ait : site_state_amplicon_stats) {
    if (ait.get_amplicon() == amplicon) {
      return &ait;
    }
  }
  site_state_amplicon_stats.emplace_back(amplicon);
  return &site_state_amplicon_stats.back();
}

void site_state_aggregator_stats::accumulate_amplicon_depths(std::unordered_map<ITNode*, uint32_t> &amp_depths) const {
  for (const auto &ait: site_state_amplicon_stats) {
    ITNode* amplicon = ait.get_amplicon();
    const site_amplicon_aggregator_stats &amplicon_stats = ait;
    amp_depths[amplicon] += amplicon_stats.get_count();
  }
}

amplicon_variation_data site_state_aggregator_stats::calculate_amplicon_variation(const std::unordered_map<ITNode *, uint32_t> &amp_depths) const {
  amplicon_variation_data result;
  std::vector<site_amplicon_aggregator_stats> unambiguous_site_state_amplicon_stats = get_unambiguous_site_state_amplicon_stats();
  result.num_amplicons = unambiguous_site_state_amplicon_stats.size();
  if (result.num_amplicons <= 1)
    return result;

  result.frequencies.reserve(result.num_amplicons);
  result.weights.reserve(result.num_amplicons);

  for(const auto& ssa_stats : unambiguous_site_state_amplicon_stats) {
    ITNode* amplicon = ssa_stats.get_amplicon();
    uint32_t count = ssa_stats.get_count();

    auto depth_it = amp_depths.find(amplicon);

    if (depth_it != amp_depths.end() && depth_it->second > 0) {
      double frequency = static_cast<double>(count) / static_cast<double>(depth_it->second);
      result.frequencies.push_back(frequency);
      result.weights.push_back(depth_it->second);
      result.amplicons.push_back(amplicon);
    }
  }

  result.stdev = calculate_weighted_standard_deviation(result.frequencies, result.weights);
  result.amplicon_frequency_variation = result.stdev > STDEV_THRESHOLD;

  return result;
}

uint32_t site_state_aggregator_stats::get_depth() const {
  uint32_t total_depth = 0;
  for (const auto& ssa_stats : site_state_amplicon_stats) {
    total_depth += ssa_stats.get_count();
  }
  return total_depth;
}

uint8_t site_state_aggregator_stats::get_mean_quality() const {
  uint8_t mean_quality;
  uint32_t weighted_sum = 0;
  uint32_t total_depth = 0;
  for (const auto& amplicon_stats : site_state_amplicon_stats) {
    uint32_t count = amplicon_stats.get_count();
    uint8_t amp_mean_quality = amplicon_stats.get_mean_quality();

    weighted_sum += count * amp_mean_quality;
    total_depth += count;
  }
  mean_quality = static_cast<uint8_t>(
      static_cast<double>(weighted_sum) / (total_depth)
  );
  return mean_quality;
}

double site_state_aggregator_stats::calculate_weighted_standard_deviation(std::vector<double> values, std::vector<uint32_t> weights) {
  double weighted_sum = 0.0, total_weight = 0.0;

  // Compute weighted mean
  for (size_t i = 0; i < values.size(); ++i) {
    weighted_sum += values[i] * weights[i];
    total_weight += weights[i];
  }
  double mean = weighted_sum / total_weight;

  // Compute weighted variance
  double variance = 0.0f;
  for (size_t i = 0; i < values.size(); ++i) {
    variance += weights[i] * std::pow(values[i] - mean, 2);
  }
  variance /= total_weight;

  return std::sqrt(variance);
}