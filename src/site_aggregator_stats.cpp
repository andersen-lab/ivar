#include "site_aggregator_stats.h"

void site_aggregator_stats::add_site(site_state ss) {
  site_state_aggregator_stats *ss_stats = find_or_create_site_state_aggregator_stats(ss.state);
  ss_stats->add_site(ss);
}

site_state_aggregator_stats* site_aggregator_stats::find_or_create_site_state_aggregator_stats(const std::string &state) {
  for(auto &site_state_stat : site_state_stats) {
    if(site_state_stat.get_state() == state) {
      return &site_state_stat;
    }
  }
  site_state_stats.emplace_back(state);
  return &site_state_stats.back();
}

const site_state_aggregator_stats* site_aggregator_stats::get_site_state_aggregator_stats(const std::string &state) const {
  for (auto &site_state_stat : site_state_stats) {
    if (site_state_stat.get_state() == state) {
      return &site_state_stat;
    }
  }
  return nullptr;
}

void site_aggregator_stats::accumulate_amplicon_depths(std::unordered_map<ITNode *, uint32_t> &amp_depths) {
  for (auto sit = site_state_stats.begin(); sit != site_state_stats.end(); ++sit) {
    if(site_state::is_insertion(sit->get_state()))
      continue;
    sit->accumulate_amplicon_depths(amp_depths);
  }
}

uint32_t site_aggregator_stats::get_depth(const std::string &state) const {
  const site_state_aggregator_stats* ss_stats = get_site_state_aggregator_stats(state);
  if(ss_stats != nullptr) {
    return ss_stats->get_depth();
  }
  return 0;
}

uint8_t site_aggregator_stats::get_mean_quality(const std::string &state) const {
  const site_state_aggregator_stats* ss_stats = get_site_state_aggregator_stats(state);
  if(ss_stats != nullptr) {
    return ss_stats->get_mean_quality();
  }
  return 0; // TODO: Add an explicit indication of null
}

double site_aggregator_stats::get_frequency(const std::string &state, uint32_t total_depth) const {
  const site_state_aggregator_stats* ss_stats = get_site_state_aggregator_stats(state);
  if (ss_stats != nullptr) {
    return ss_stats->get_depth() / static_cast<double>(total_depth);
  }
  return 0;
}

// Ungapped depth
uint32_t site_aggregator_stats::get_total_depth() const {
  uint32_t total_depth = 0;
  for (const auto& sit : site_state_stats) {
    std::string state = sit.get_state();
    if(site_state::is_insertion(state) || site_state::is_deletion(state) || site_state::is_gap(state))
      continue;
    total_depth += sit.get_depth();
  }
  return total_depth;
}

double site_aggregator_stats::get_gapped_frequency(const std::string &state, uint32_t total_gapped_depth) const {
  const site_state_aggregator_stats* ss_stats = get_site_state_aggregator_stats(state);
  if (ss_stats != nullptr) {
    return ss_stats->get_depth() / static_cast<double>(total_gapped_depth);
  }
  return 0;
}

// Gapped depth
uint32_t site_aggregator_stats::get_total_gapped_depth() const {
  uint32_t total_gapped_depth = 0;
  for (const auto& sit : site_state_stats) {
    std::string state = sit.get_state();
    if(site_state::is_insertion(state))
      continue;
    total_gapped_depth += sit.get_depth();
  }
  return total_gapped_depth;
}


