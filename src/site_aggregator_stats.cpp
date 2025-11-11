#include "site_aggregator_stats.h"

void site_aggregator_stats::add_site(site_state ss) {
  site_amplicon_stats[ss.amplicon].add_site(ss);
}

void site_aggregator_stats::calculate_amplicon_stats() {

}