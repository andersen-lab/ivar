#include "site_aggregator_stats.h"

void site_aggregator_stats::add_site(site_state ss) {
  site_state_stats[ss.state].add_site(ss);
}

void site_aggregator_stats::accumulate_amplicon_depths(std::unordered_map<ITNode *, uint32_t> &amp_depths) {
  for (auto sit = site_state_stats.begin(); sit != site_state_stats.end(); ++sit) {
    const std::string &state = sit->first;
    const site_state_aggregator_stats &state_stats = sit->second;
    if(site_state::is_insertion(state))
      continue;
    state_stats.accumulate_amplicon_depths(amp_depths);
  }
}

uint32_t site_aggregator_stats::get_depth(const std::string &state) const {
  auto it = site_state_stats.find(state);
  if (it != site_state_stats.end()) {
    return it->second.get_depth();
  }
  return 0;
}

uint8_t site_aggregator_stats::get_mean_quality(const std::string &state) const {
  auto it = site_state_stats.find(state);
  if (it != site_state_stats.end()) {
    return it->second.get_mean_quality();
  }
  return 0; // TODO: Add an explicit indication of null
}

double site_aggregator_stats::get_frequency(const std::string &state, uint32_t total_depth) const {
  auto it = site_state_stats.find(state);
  if (it != site_state_stats.end()) {
    return it->second.get_depth() / static_cast<double>(total_depth);
  }
  return 0;
}

// Ungapped depth
uint32_t site_aggregator_stats::get_total_depth() const {
  uint32_t total_depth = 0;
  for (const auto& sit : site_state_stats) {
    if(site_state::is_insertion(sit.first) || site_state::is_deletion(sit.first) || site_state::is_gap(sit.first))
      continue;
    total_depth += sit.second.get_depth();
  }
  return total_depth;
}

double site_aggregator_stats::get_gapped_frequency(const std::string &state, uint32_t total_gapped_depth) const {
  auto it = site_state_stats.find(state);
  if (it != site_state_stats.end()) {
    return it->second.get_depth() / static_cast<double>(total_gapped_depth);
  }
  return 0;
}

// Gapped depth
uint32_t site_aggregator_stats::get_total_gapped_depth() const {
  uint32_t total_gapped_depth = 0;
  for (const auto& sit : site_state_stats) {
    if(site_state::is_insertion(sit.first))
      continue;
    total_gapped_depth += sit.second.get_depth();
  }
  return total_gapped_depth;
}


