#include "site_aggregator.h"

const std::unordered_map<site_aggregator_key, site_aggregator_stats>& site_aggregator::get_data() {
  return aggregated_site_states;
}

void site_aggregator::aggregate(const std::vector<site_state> site_states) {
  for(int i = 0;i < site_states.size(); i++){
    // Aggregate site level stats by position, state, and amplicon
    site_aggregator_key key = {
      site_states[i].coordinate,
      site_states[i].state
    };
    aggregated_site_states[key].add_site(site_states[i]);
  }
}
