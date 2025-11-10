#ifndef HASH_UTILS_H
#define HASH_UTILS_H

#include <cstddef>

// From https://stackoverflow.com/questions/4948780/magic-number-in-boosthash-combine
namespace hash_utils {
  inline size_t hash_combine(size_t seed, size_t value) {
    return seed ^ (value + 0x9e3779b9 + (seed << 6) + (seed >> 2));
  }
}

#endif