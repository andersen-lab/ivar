#include "site_state.h"

const std::string site_state::GAP = "_";

void site_state::set_nucleotide(const std::string &nucleotide, uint8_t qual, uint32_t position) {
  state = nucleotide;
  quality = qual;
  coordinate.type = NUCLEOTIDE;
  coordinate.position = position;
}

void site_state::set_nucleotide(const char &nucleotide, uint8_t qual, uint32_t position) {
  state.assign(1, nucleotide);
  quality = qual;
  coordinate.type = NUCLEOTIDE;
  coordinate.position = position;
}

void site_state::set_amplicon(ITNode* node, bool is_ambiguous) {
  amplicon = node;
  this->is_ambiguous = is_ambiguous;
}

void site_state::set_nucleotide_gap(uint8_t min_qual, uint32_t position) {
  state = GAP;
  quality = min_qual;
  coordinate.type = NUCLEOTIDE;
  coordinate.position = position;
}
bool site_state::is_deletion(const std::string &state) {
    return state[0] == '-';
}

bool site_state::is_insertion(const std::string &state) {
    return state[0] == '+';
}

bool site_state::is_gap(const std::string &state) {
    return state == GAP;
}