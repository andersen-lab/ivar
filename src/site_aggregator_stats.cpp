//
// Created by Karthik on 11/9/25.
//

#include "site_aggregator_stats.h"

void site_aggregator_stats::add_site(site_state ss) {
  bool is_del = false;
  bool is_ins = false;

  if (ss.state[0] == '-') {
    is_del = true;
  } else if (ss.state[0] == '+') {
    is_ins = true;
  }

  if (is_del) {
    gapped_depth += 1;
    mean_quality = static_cast<uint8_t>(
        mean_quality + static_cast<double>(ss.quality - mean_quality) / gapped_depth
    );
  } else if (is_ins) {
    depth += 1;
    mean_quality = static_cast<uint8_t>(
        mean_quality + static_cast<double>(ss.quality - mean_quality) / gapped_depth
    );
    gapped_depth += 1;
  } else {
    depth += 1;
    gapped_depth += 1;
    mean_quality = static_cast<uint8_t>(
        mean_quality + static_cast<double>(ss.quality - mean_quality) / depth
    );
  }
}