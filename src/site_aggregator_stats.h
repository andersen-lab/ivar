#include <unordered_map>
#include "interval_tree.h"
#include "site_amplicon_aggregator_stats.h"

#ifndef SITE_AGGREGATOR_STATS_H
#define SITE_AGGREGATOR_STATS_H

class site_aggregator_stats {
 private:
  std::unordered_map<ITNode*, site_amplicon_aggregator_stats> site_amplicon_stats;
 public:
  void add_site(site_state ss);
  void calculate_amplicon_stats();
  std::unordered_map<ITNode*, site_amplicon_aggregator_stats> get_site_amplicon_stats() {
    return site_amplicon_stats;
  }
};

#endif  // SITE_AGGREGATOR_STATS_H
