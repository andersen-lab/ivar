#include "site_aggregator.h"

const std::vector<site_aggregator_stats>& site_aggregator::get_data() {
  return aggregated_site_states;
}

void site_aggregator::aggregate(const std::vector<site_state> &site_states) {
  for(int i = 0;i < site_states.size(); i++){
    // Aggregate site level stats by position
    aggregated_site_states[site_states[i].coordinate.position].add_site(site_states[i]);
  }
}

bool site_aggregator::calculate_amplicon_depths(site_coordinate coord, std::unordered_map<ITNode*, uint32_t> &amp_depths) {
  amp_depths.clear();
  if(aggregated_site_states.size() < coord.position)
    return false;
  aggregated_site_states[coord.position].accumulate_amplicon_depths(amp_depths);
  return true;
}

void site_aggregator::clear() {
  aggregated_site_states.clear();
}
