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

  site_aggregator sa;
  sa.initialize(24);
  sa.aggregate(site_states);
  std::vector<site_aggregator_stats> aggregated_site_states = sa.get_data();

  // Expected values
  std::vector<site_coordinate> expected_coord_keys = {
    {
      NUCLEOTIDE,
      1
    },
    {
      NUCLEOTIDE,
      1
    },
    {
      NUCLEOTIDE,
      22
    },
    {
      NUCLEOTIDE,
      23
    },
    {
      NUCLEOTIDE,
      23
    },
    {
      NUCLEOTIDE,
      24
    }
  };

  std::vector<std::string> expected_site_states = {
    "A",
    "T",
    "-TC",
    site_state::GAP,
    "C",
    "+G"
  };

  std::vector<site_amplicon_aggregator_stats> expected_site_amplicon_aggregator_stats;
  site_amplicon_aggregator_stats saas(amp1[0]);
  // 1 A
  saas.set_stats(1, 30);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 1 T
  saas.clear();
  saas.set_stats(2, 35);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 22 -TC
  saas.clear();
  saas.set_stats(2, 22);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 23 GAP
  saas.clear();
  saas.set_stats(2, 20);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 23 C
  saas.clear();
  saas.set_stats(1, 20);
  expected_site_amplicon_aggregator_stats.push_back(saas);
  // 24 +G
  saas.clear();
  saas.set_stats(2, 20);
  expected_site_amplicon_aggregator_stats.push_back(saas);

  int pass = 0;
  for(int i = 0; i < expected_site_amplicon_aggregator_stats.size(); i++){
    if(*(aggregated_site_states[expected_coord_keys[i].position]
            .get_site_state_aggregator_stats(expected_site_states[i])
            ->get_site_amplicon_aggregator_stats(amp1[0])) == expected_site_amplicon_aggregator_stats[i]) {
      pass +=1;
    } else {
      std::cerr << "Aggregated stats not match for site " << expected_coord_keys[i].position << " " << expected_site_states[i] << std::endl;
      std::cerr << "Observed states for site " << expected_coord_keys[i].position << ":" << aggregated_site_states[expected_coord_keys[i].position].get_site_state_aggregator_stats(expected_site_states[i])->get_site_amplicon_aggregator_stats(amp1[0])->get_count() << std::endl;
    }
  }

  if(pass != expected_site_amplicon_aggregator_stats.size()) {
    std::cerr << "Failed site amplicon aggregator stats test." << std::endl;
    return 1;
  }

  site_states.clear();

  // Position 1
  set_site_state("A", 30, 1, amp2[0], false, site_states);
  set_site_state("T", 25, 1, amp1[0], false, site_states);
  set_site_state("T", 45, 1, amp1[0], false, site_states);

  // Position 22
  set_site_state("-TC", 22, 22, amp1[0], false, site_states);
  set_site_state("-TC", 22, 22, amp1[0], false, site_states);

  // Position 23
  // Gap from -TC at 22
  set_site_state_nucleotide_gap(20, 23, amp1[0], false, site_states);
  set_site_state_nucleotide_gap(20, 23, amp2[0], false, site_states);

  // Position 24
  set_site_state("+G", 20, 24, amp1[0], false, site_states);
  set_site_state("+G", 20, 24, amp2[0], false, site_states);

  sa.clear();
  aggregated_site_states.clear();

  sa.initialize(24);
  sa.aggregate(site_states);

  // Expected amplicon depths
  std::vector<std::unordered_map<ITNode*, uint32_t>> expected_amp_depths = {
      {
        {amp1[0], 2},
        {amp2[0], 1}
      },
      {
        {amp1[0], 2}
      },
      {
        {amp1[0], 1},
        {amp2[0], 1}
      },
      {}
  };

  std::unordered_map<ITNode*, uint32_t> amp_depths;
  site_coordinate coord{
  NUCLEOTIDE,
  1
  };

  int positions[] = {1, 22, 23, 24};

  for(int i = 0; i < 4; i++) {
    coord.position = positions[i];
    amp_depths.clear();
    sa.calculate_amplicon_depths(coord, amp_depths);
    if(amp_depths.size() != expected_amp_depths[i].size()) {
      for(auto it = amp_depths.begin(); it != amp_depths.end(); ++it){
        std::cerr << it->second << std::endl;
      }
      std::cerr << "Amplicon depths size mismatch for " << i << std::endl;
      return 1;
    }

    for (auto it = expected_amp_depths[i].begin(); it != expected_amp_depths[i].end(); ++it){
      ITNode *amplicon = it->first;
      uint32_t depth = it->second;
      if(amp_depths[amplicon] != depth) {
        std::cerr << "Amplicon depth mismatch for amplicon: [" << amplicon->data->low << ", " << amplicon->data->high << "]" << std::endl;
        return 1;
      }
    }


  }
  return 0;
}