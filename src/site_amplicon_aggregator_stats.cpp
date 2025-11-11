#include "site_amplicon_aggregator_stats.h"

void site_amplicon_aggregator_stats::add_site(site_state ss) {
  //Note: site_amplicon_aggregator_stats stores for one state and one amplicon
  //TODO: Add an explicit check that ss.state and ss.amplicon match the stored ones
  if (ss.state[0] == '-' || ss.state == site_state::GAP) {
    // Deletion
    gapped_depth += 1;
    mean_quality = static_cast<uint8_t>(
        mean_quality + static_cast<double>(ss.quality - mean_quality) / gapped_depth
    );
  } else if (ss.state[0] == '+') {
    // Insertion
    depth += 1;
    gapped_depth += 1;
    mean_quality = static_cast<uint8_t>(
        mean_quality + static_cast<double>(ss.quality - mean_quality) / gapped_depth
    );
  } else {
    // SNVs
    depth += 1;
    gapped_depth += 1;
    mean_quality = static_cast<uint8_t>(
        mean_quality + static_cast<double>(ss.quality - mean_quality) / depth
    );
  }
}
bool site_amplicon_aggregator_stats::clear() {
  gapped_depth = 0;
  depth = 0;
  mean_quality = 0;
  return true;
}
