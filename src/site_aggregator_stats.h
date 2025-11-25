#include <unordered_map>
#include "interval_tree.h"
#include "site_amplicon_aggregator_stats.h"
#include "site_state_aggregator_stats.h"

#ifndef SITE_AGGREGATOR_STATS_H
#define SITE_AGGREGATOR_STATS_H

class site_aggregator_stats {
 private:
  // Map from state to its stats
  std::unordered_map<std::string, site_state_aggregator_stats> site_state_stats;
 public:
  void add_site(site_state ss);
  std::unordered_map<std::string, site_state_aggregator_stats> get_site_state_stats() const {
    return site_state_stats;
  }
  void accumulate_amplicon_depths(std::unordered_map<ITNode*, uint32_t> &amp_depths);

  // Get depth for a given state
  uint32_t get_depth(const std::string &state) const;

  // Get mean quality for a given state
  uint8_t get_mean_quality(const std::string &state) const;

  // Get frequency for a given state
  double get_frequency(const std::string &state, uint32_t total_depth) const;

  // Get total depth across all states
  uint32_t get_total_depth() const;

  // Get gapped frequency for a given state
  double get_gapped_frequency(const std::string &state, uint32_t total_gapped_depth) const;

  // Get total gapped depth across all states at a coordinate
  uint32_t get_total_gapped_depth() const;
};

#endif  // SITE_AGGREGATOR_STATS_H
