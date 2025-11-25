#include <cstdint>
#include "site_state.h"

#ifndef SITE_AMPLICON_AGGREGATOR_STATS_H
#define SITE_AMPLICON_AGGREGATOR_STATS_H

class site_amplicon_aggregator_stats {
 private:
  uint32_t count=0;
  uint8_t mean_quality=0;

 public:
  void add_site(site_state ss);

  void set_stats(uint32_t c, uint8_t mq) {
    count = c;
    mean_quality = mq;
  }

  uint32_t get_count() const {
    return count;
  }

  uint8_t get_mean_quality() const {
    return mean_quality;
  }

  bool operator == (const site_amplicon_aggregator_stats &other) const {
    return (count == other.count);
  }

  bool clear();
};

#endif  // SITE_AMPLICON_AGGREGATOR_STATS_H
