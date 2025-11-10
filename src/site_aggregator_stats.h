#include <cstdint>
#include "site_state.h"

#ifndef SITE_AGGREGATOR_STATS_H
#define SITE_AGGREGATOR_STATS_H

class site_aggregator_stats {
 private:
  uint32_t pos;
  uint32_t depth=0;
  uint32_t gapped_depth=0;
  uint32_t mean_quality=0;

 public:
  void add_site(site_state ss);
  uint32_t get_pos() const {return pos;};
  uint32_t get_depth() const { return depth; };
  uint32_t get_gapped_depth() const { return gapped_depth; };
  uint32_t get_qual() const { return mean_quality; };
};

#endif  // SITE_AGGREGATOR_STATS_H
