#include "../src/site_amplicon_aggregator_stats.h"

void set_site_state(std::string nuc, uint8_t quality, uint32_t position, ITNode *amp, bool is_ambiguous, std::vector<site_state> &expected_states) {
  site_state ss;
  ss.set_nucleotide(nuc, quality, position);
  ss.set_amplicon(amp, is_ambiguous);
  expected_states.push_back(ss);
}

void set_site_amplicon_aggregator_stats(std::vector<site_state> &site_states, site_amplicon_aggregator_stats &saas) {
  for(int i = 0; i < site_states.size(); i++){
    saas.add_site(site_states[i]);
  }
}

site_amplicon_aggregator_stats get_expected_site_amplicon_aggregator_stats(uint32_t depth, uint8_t mean_quality) {
  site_amplicon_aggregator_stats expected_sas;
  expected_sas.set_stats(depth, mean_quality);
  return expected_sas;
}

void set_site_state_nucleotide_gap(uint8_t quality, uint32_t position, ITNode *amp, bool is_ambiguous, std::vector<site_state> &expected_states) {
  site_state ss;
  ss.set_nucleotide_gap(quality, position);
  ss.set_amplicon(amp, is_ambiguous);
  expected_states.push_back(ss);
}

int main() {
  std::vector<site_state> site_states;

  int pass = 0;

  // Position 1 state A
  set_site_state("A", 30, 1, nullptr, false, site_states);
  site_amplicon_aggregator_stats saas;
  set_site_amplicon_aggregator_stats(site_states, saas);
  if (saas == get_expected_site_amplicon_aggregator_stats(1, 30)){
      pass += 1;
  } else {
    std::cerr << "Failed at position 1 state A" << std::endl;
    return 1;
  }

  // Position 1 state T
  site_states.clear();
  saas.clear();
  set_site_state("T", 25, 1, nullptr, false, site_states);
  set_site_state("T", 45, 1, nullptr, false, site_states);

  set_site_amplicon_aggregator_stats(site_states, saas);

  if (saas == get_expected_site_amplicon_aggregator_stats(2, 35)){
    pass += 1;
  } else {
    std::cerr << "Failed at position 1 state T" << std::endl;
    return 1;
  }

  // Position 22 state -TC
  site_states.clear();
  saas.clear();
  set_site_state("-TC", 22, 22, nullptr, false, site_states);
  set_site_state("-TC", 24, 22, nullptr, false, site_states);

  set_site_amplicon_aggregator_stats(site_states, saas);

  if (saas == get_expected_site_amplicon_aggregator_stats(2, 23)){
    pass += 1;
  } else {
    std::cerr << "Failed at position 22 state -TC" << std::endl;
    return 1;
  }

  // Position 23 GAP
  // Gap from -TC at 22
  site_states.clear();
  saas.clear();
  set_site_state_nucleotide_gap(20, 23, nullptr, false, site_states);
  set_site_state_nucleotide_gap(21, 23, nullptr, false, site_states);

  set_site_amplicon_aggregator_stats(site_states, saas);

  if (saas == get_expected_site_amplicon_aggregator_stats(2, 20)){
    pass += 1;
  } else {
    std::cerr << "Failed at position 23 state GAP" << std::endl;
    return 1;
  }

  // Position 23 state C
  site_states.clear();
  saas.clear();
  set_site_state("C", 20, 23, nullptr, false, site_states);

  set_site_amplicon_aggregator_stats(site_states, saas);

  if (saas == get_expected_site_amplicon_aggregator_stats(1, 20)){
    pass += 1;
  } else {
    std::cerr << "Failed at position 23 state C" << std::endl;
    return 1;
  }

  // Position 24 state +G
  site_states.clear();
  saas.clear();

  set_site_state("+G", 20, 24, nullptr, false, site_states);
  set_site_state("+G", 20, 24, nullptr, false, site_states);

  set_site_amplicon_aggregator_stats(site_states, saas);

  if (saas == get_expected_site_amplicon_aggregator_stats(2, 20)){
    pass += 1;
  } else {
    std::cerr << "Failed at position 24 state +G" << std::endl;
    return 1;
  }

  return (pass == 6) ? 0 : 1;
}