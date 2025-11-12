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


  site_aggregator sa;
  sa.aggregate(site_states);
  std::unordered_map<site_aggregator_key, site_aggregator_stats> aggregated_site_states = sa.get_data();

  std::vector<site_aggregator_key> expected_keys = {
    {
      {
          NUCLEOTIDE,
          1
      },
      "A"
    },
    {
      {
          NUCLEOTIDE,
          1
      },
      "T"
    },
    {
      {
          NUCLEOTIDE,
          22
      },
      "-TC"
    },
    {
      {
      NUCLEOTIDE,
      23
      },
      site_state::GAP
    },
    {
      {
      NUCLEOTIDE,
      23
      },
      "C"
    },
    {
      {
      NUCLEOTIDE,
      24
      },
      "+G"
    }
  };

  std::vector<site_amplicon_aggregator_stats> expected_site_amplicon_aggregator_stats;
  site_amplicon_aggregator_stats saas;
  // 1 A
  saas.set_stats(1, 1, 30);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 1 T
  saas.clear();
  saas.set_stats(2, 2, 35);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 22 -TC
  saas.clear();
  saas.set_stats(0, 2, 22);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 23 GAP
  saas.clear();
  saas.set_stats(0, 2, 20);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 23 C
  saas.clear();
  saas.set_stats(1, 1, 20);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 24 +G
  saas.clear();
  saas.set_stats(2, 2, 20);
  expected_site_amplicon_aggregator_stats.push_back(saas);

  int pass = 0;
  for(int i = 0; i < expected_site_amplicon_aggregator_stats.size(); i++){
    if(aggregated_site_states[expected_keys[i]].get_site_amplicon_stats()[amp1[0]] == expected_site_amplicon_aggregator_stats[i]) {
      pass +=1;
    } else {
      std::cerr << "Aggregated stats not match for site " << expected_keys[i].coordinate.position << " " << expected_keys[i].state << std::endl;
    }
  }

  return (pass == expected_site_amplicon_aggregator_stats.size()) ? 0 : 1;
}