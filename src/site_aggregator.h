#include <cstdint>
#include <unordered_map>
#include "site_state.h"
#include "site_aggregator_stats.h"
#include <functional>
#include "hash_utils.h"

#ifndef SITE_AGGREGATOR_H
#define SITE_AGGREGATOR_H

struct site_aggregator_key {
  site_coordinate coordinate;
  std::string state;
  ITNode* amplicon;

  size_t hash() const {
    size_t h = coordinate.hash();
    h = hash_utils::hash_combine(h, std::hash<std::string>{}(state));
    if (amplicon != nullptr) {
      h = hash_utils::hash_combine(h, static_cast<size_t>(amplicon->data->low));
      h = hash_utils::hash_combine(h, static_cast<size_t>(amplicon->data->high));
    }
    return h;
  }

  bool operator==(const site_aggregator_key& other) const {
    if (!(coordinate == other.coordinate && state == other.state)) {
      return false;
    }

    if (amplicon == nullptr || other.amplicon == nullptr) {
      return amplicon == other.amplicon;
    }

    return amplicon->data->low == other.amplicon->data->low &&
           amplicon->data->high == other.amplicon->data->high;
  }

};

namespace std {
template<>
struct hash<site_aggregator_key> {
  size_t operator()(const site_aggregator_key& c) const {
    return c.hash();
  }
};
}

class site_aggregator {
 private:
  std::unordered_map<site_aggregator_key, site_aggregator_stats> aggregated_site_states;
 public:
  void aggregate(const std::vector<site_state> site_states);
  const std::unordered_map<site_aggregator_key, site_aggregator_stats>& get_data();
};

#endif  // SITE_AGGREGATOR_H
