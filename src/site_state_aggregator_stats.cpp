#include "site_state_aggregator_stats.h"
#include <algorithm>
#include <cmath>

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

  // Every amplicon covering the position, not just those carrying this state, so
  // that a state absent from an amplicon is reported at frequency 0 rather than
  // dropped
  for (const auto &depth_it : amp_depths) {
    ITNode* amplicon = depth_it.first;
    uint32_t depth = depth_it.second;
    if (amplicon == nullptr || depth == 0)
      continue;
    result.amplicons.push_back(amplicon);
  }

  std::sort(result.amplicons.begin(), result.amplicons.end(), [](ITNode* a, ITNode* b) {
    if (a->data->low != b->data->low)
      return a->data->low < b->data->low;
    return a->data->high < b->data->high;
  });

  result.num_amplicons = result.amplicons.size();
  result.frequencies.reserve(result.num_amplicons);
  result.weights.reserve(result.num_amplicons);

  for (ITNode* amplicon : result.amplicons) {
    uint32_t depth = amp_depths.at(amplicon);
    const site_amplicon_aggregator_stats* ssa_stats = get_site_amplicon_aggregator_stats(amplicon);
    uint32_t count = (ssa_stats != nullptr) ? ssa_stats->get_count() : 0;

    result.frequencies.push_back(static_cast<double>(count) / static_cast<double>(depth));
    result.weights.push_back(depth);
  }

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
  double weighted_sum = 0;
  uint32_t total_depth = 0;
  for (const auto& amplicon_stats : site_state_amplicon_stats) {
    uint32_t count = amplicon_stats.get_count();
    double amp_mean_quality = amplicon_stats.get_mean_quality();

    weighted_sum += count * amp_mean_quality;
    total_depth += count;
  }
  mean_quality = static_cast<uint8_t>(
      weighted_sum / (total_depth)
  );
  return mean_quality;
}