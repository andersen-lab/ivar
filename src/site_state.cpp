#include "site_state.h"

void site_state::set_nucleotide(std::string nucleotide, uint8_t qual, uint32_t position) {
  state = nucleotide;
  quality = qual;
  coordinate.type = NUCLEOTIDE;
  coordinate.position = position;
}

void site_state::set_amplicon(ITNode* node, bool is_ambiguous) {
  amplicon = node;
  this->is_ambiguous = is_ambiguous;
}