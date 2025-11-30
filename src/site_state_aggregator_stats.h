#include "site_amplicon_aggregator_stats.h"
#include "site_state.h"
#include <unordered_map>

#ifndef SITE_STATE_AGGREGATOR_STATS_H
#define SITE_STATE_AGGREGATOR_STATS_H

struct amplicon_variation_data {
  std::vector<double> frequencies;
  std::vector<uint32_t> weights;
  std::vector<ITNode*> amplicons;

  int num_amplicons = 0;
  double stdev = 0;
  bool amplicon_frequency_variation = false;

  bool has_variation() const {
    return num_amplicons > 1;
  }
};

class site_state_aggregator_stats {
 private:
  // Map from base at position to its amplicon stats
  std::vector<site_amplicon_aggregator_stats> site_state_amplicon_stats;
  std::string state;

  std::vector<site_amplicon_aggregator_stats> get_unambiguous_site_state_amplicon_stats() const {
    std::vector<site_amplicon_aggregator_stats> filtered;
    for (const auto& ssa_stats : site_state_amplicon_stats) {
      if (ssa_stats.get_amplicon() != nullptr) {
        filtered.emplace_back(ssa_stats);
      }
    }
    return filtered;
  }

  static double calculate_weighted_standard_deviation(std::vector<double> values, std::vector<uint32_t> weights);
  site_amplicon_aggregator_stats* find_or_create_site_amplicon_aggregator_stats(ITNode* amplicon);

  static const double STDEV_THRESHOLD;

 public:

  const std::string& get_state() const {
    return state;
  }

  const site_amplicon_aggregator_stats* get_site_amplicon_aggregator_stats(ITNode* amplicon) const;

  explicit site_state_aggregator_stats(const std::string &state){
    this->state = state;
  }

  void add_site(site_state &ss);

  std::vector<site_amplicon_aggregator_stats> get_amplicon_stats() const {
    return site_state_amplicon_stats;
  }

  void accumulate_amplicon_depths(std::unordered_map<ITNode*, uint32_t>& amp_depths) const;
  amplicon_variation_data calculate_amplicon_variation(const std::unordered_map<ITNode*, uint32_t>& amp_depths) const;

  uint32_t get_depth() const;
  uint8_t get_mean_quality() const;


};

#endif  // SITE_STATE_AGGREGATOR_STATS_H
