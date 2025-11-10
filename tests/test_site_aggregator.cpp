#include "../src/site_aggregator.h"
#include <unordered_map>

void set_site_state(std::string nuc, uint8_t quality, uint32_t position, ITNode *amp, bool is_ambiguous, std::vector<site_state> &expected_states) {
  site_state ss;
  ss.set_nucleotide(nuc, quality, position);
  ss.set_amplicon(amp, is_ambiguous);
  expected_states.push_back(ss);
}

void set_site_state_nucleotide_gap(uint8_t quality, uint32_t position, ITNode *amp, bool is_ambiguous, std::vector<site_state> &expected_states) {
  site_state ss;
  ss.set_nucleotide_gap(quality, position);
  ss.set_amplicon(amp, is_ambiguous);
  expected_states.push_back(ss);
}

void set_expected_stats(uint32_t position, std::string nuc, uint32_t depth, uint32_t gapped_depth, uint8_t mean_quality, ITNode *amp, std::unordered_map<site_aggregator_key, site_aggregator_stats> &expected_aggregated_site_states) {
  site_aggregator_stats expected_sas;
  site_aggregator_key expected_key;
  expected_sas.set_stats(depth, gapped_depth, mean_quality);
  expected_key = {
      {
          NUCLEOTIDE,
          position
      },
      nuc,
      amp
  };
  expected_aggregated_site_states[expected_key] = expected_sas;
}

int main() {
  IntervalTree test_amplicons = IntervalTree();

  test_amplicons.insert(Interval(1, 5));
  test_amplicons.insert(Interval(5, 30));

  std::vector<ITNode*> amp1;
  test_amplicons.find_read_amplicon(1, 2, amp1);

  std::vector<ITNode*> amp2;
  test_amplicons.find_read_amplicon(20, 25, amp2);

  std::vector<site_state> site_states;
  site_state ss;

  // Position 1
  set_site_state("A", 30, 1, amp1[0], false, site_states);
  set_site_state("T", 25, 1, amp1[0], false, site_states);
  set_site_state("T", 45, 1, amp1[0], false, site_states);

  // Position 22
  set_site_state("-TC", 22, 22, amp1[0], false, site_states);
  set_site_state("-TC", 22, 22, amp1[0], false, site_states);

  // Position 23
  // Gap from -TC at 22
  set_site_state_nucleotide_gap(20, 23, amp1[0], false, site_states);
  set_site_state_nucleotide_gap(20, 23, amp1[0], false, site_states);
  set_site_state("C", 20, 23, amp1[0], false, site_states);

  // Position 24
  set_site_state("+G", 20, 24, amp1[0], false, site_states);
  set_site_state("+G", 20, 24, amp1[0], false, site_states);

  // Expected values
  std::unordered_map<site_aggregator_key, site_aggregator_stats> expected_aggregated_site_states;

  set_expected_stats(1, "A", 1, 1, 30, amp1[0], expected_aggregated_site_states);
  set_expected_stats(1, "T", 2, 2, 35, amp1[0], expected_aggregated_site_states);

  set_expected_stats(22, "-TC", 0, 2, 22, amp1[0], expected_aggregated_site_states);

  set_expected_stats(23, site_state::GAP, 0, 2, 20, amp1[0], expected_aggregated_site_states);
  set_expected_stats(23, "C", 1, 1, 20, amp1[0], expected_aggregated_site_states);

  set_expected_stats(24, "+G", 2, 2, 20, amp1[0], expected_aggregated_site_states);

  site_aggregator sa;
  sa.aggregate(site_states);
  std::unordered_map<site_aggregator_key, site_aggregator_stats> aggregated_site_states = sa.get_data();
  int pass = 0;

//  site_aggregator_key expected_key = {
//      {
//          NUCLEOTIDE,
//          22
//      },
//      "-TC",
//      amp1[0]
//  };
//
//  std::cout << aggregated_site_states[expected_key].get_depth() << " "
//            << aggregated_site_states[expected_key].get_gapped_depth() << " "
//            << static_cast<int>(aggregated_site_states[expected_key].get_qual()) << std::endl;
//
//  std::cout << expected_aggregated_site_states[expected_key].get_depth() << " "
//            << expected_aggregated_site_states[expected_key].get_gapped_depth() << " "
//            << static_cast<int>(expected_aggregated_site_states[expected_key].get_qual()) << std::endl;

  for (const auto& [key, value] : aggregated_site_states) {
    if(expected_aggregated_site_states[key] == value) {
      pass += 1;
    } else {
      std::cerr << "Aggregated stats not match for site " << key.coordinate.position << " " << key.state;
      if(key.amplicon != nullptr)
        std::cerr << " (" << key.amplicon->data->low << ", " << key.amplicon->data->high << ")" << std::endl;
      else
        std::cerr << " (, )" << std::endl;
      std::cerr << value.get_depth() << " " << value.get_gapped_depth() << " " << static_cast<int>(value.get_qual()) << std::endl << std::endl;
    }
  }

  if(pass != expected_aggregated_site_states.size())
    std::cerr << "Site aggregation did not work for all sites" << std::endl;
  return 0;
}