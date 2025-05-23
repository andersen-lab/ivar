#include <string>
#include <vector>

int check_position_exists(uint32_t p, std::vector<position> positions) {
  for (uint32_t i=0; i < positions.size(); i++) {
    if (p == positions[i].pos) {
      return((int)i);
    }
  }
  return(-1);
}
