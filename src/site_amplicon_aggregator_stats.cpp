#include "site_amplicon_aggregator_stats.h"

void site_amplicon_aggregator_stats::add_site(uint8_t ss_quality) {
  //TODO: Add explicit check if site_state.amplicon == this->amplicon
  count += 1;
  mean_quality = static_cast<uint8_t>(
      mean_quality + static_cast<double>(ss_quality - mean_quality) / count
  );
}

bool site_amplicon_aggregator_stats::clear() {
  count = 0;
  return true;
}
