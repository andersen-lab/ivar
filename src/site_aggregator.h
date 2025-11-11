#include <cstdint>
#include <functional>
#include <unordered_map>

#include "hash_utils.h"
#include "site_aggregator_stats.h"
#include "site_state.h"

#ifndef SITE_AGGREGATOR_H
#define SITE_AGGREGATOR_H

struct site_aggregator_key {
  site_coordinate coordinate;
  std::string state;

  size_t hash() const {
    size_t h = coordinate.hash();
    h = hash_utils::hash_combine(h, std::hash<std::string>{}(state));
    return h;
  }

  bool operator==(const site_aggregator_key& other) const {
    return (coordinate == other.coordinate) && (state == other.state);
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
