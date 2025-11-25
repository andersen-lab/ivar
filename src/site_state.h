#include <cstdint>
#include <string>
#include "interval_tree.h"
#include "hash_utils.h"

#ifndef SITE_STATE_H
#define SITE_STATE_H

enum site_type {
  NUCLEOTIDE,
  AMINO_ACID
};

struct site_coordinate {
  site_type type;
  uint32_t position;

  size_t hash() const {
    size_t h = static_cast<size_t>(type);
    h = hash_utils::hash_combine(h, static_cast<size_t>(position));
    return h;
  }

  bool operator == (const site_coordinate &other) const {
    return (type == other.type) && (position == other.position);
  }

  bool operator > (const site_coordinate &other) const {
    if (type != other.type) throw std::logic_error("cross-type comparison");
    return position > other.position;
  }

  bool operator < (const site_coordinate &other) const {
    if (type != other.type) throw std::logic_error("cross-type comparison");
    return position < other.position;
  }

  bool operator >= (const site_coordinate &other) const {
    if (type != other.type) throw std::logic_error("cross-type comparison");
    return position >= other.position;
  }

  bool operator <= (const site_coordinate &other) const {
    if (type != other.type) throw std::logic_error("cross-type comparison");
    return position <= other.position;
  }
};

namespace std {
  template<>
  struct hash<site_coordinate> {
    size_t operator()(const site_coordinate& c) const {
      return c.hash();
    }
  };
}

class site_state {
 public:
  site_coordinate coordinate;
  std::string state;
  uint8_t quality;
  static const std::string GAP;
  static bool is_deletion(const std::string &state);
  static bool is_insertion(const std::string &state);
  static bool is_gap(const std::string &state);

  ITNode* amplicon;
  bool is_ambiguous;

  void set_nucleotide(std::string nucleotide, uint8_t qual, uint32_t position);
  void set_amplicon(ITNode* node, bool is_ambiguous = false);
  void set_nucleotide_gap(uint8_t min_qual, uint32_t position);

  bool operator == (const site_state &other) const {
    return (coordinate == other.coordinate) && (state == other.state) && (quality == other.quality);
  }
};

#endif  // SITE_STATE_H
