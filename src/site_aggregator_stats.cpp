//
// Created by Karthik on 11/9/25.
//

#include "site_aggregator_stats.h"

void site_aggregator_stats::add_site(site_state ss) {
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