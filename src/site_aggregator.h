#include <cstdint>
#include <functional>
#include <unordered_map>

#include "hash_utils.h"
#include "site_aggregator_stats.h"
#include "site_state.h"

#ifndef SITE_AGGREGATOR_H
#define SITE_AGGREGATOR_H

class site_aggregator {
 private:
  std::vector<site_aggregator_stats> aggregated_site_states;

 public:
  bool initialize(int64_t ref_len) {
    if(ref_len >= 0) {
        aggregated_site_states.resize(ref_len + 1);  // To account for 1-based position
        return true;
    }
    return false;
  }
  void aggregate(const std::vector<site_state> &site_states);
  const std::vector<site_aggregator_stats>& get_data();
  bool calculate_amplicon_depths(site_coordinate coord, std::unordered_map<ITNode*, uint32_t> &amp_depths);
  void clear();
};

#endif  // SITE_AGGREGATOR_H
