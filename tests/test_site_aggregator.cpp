#include "../src/site_aggregator.h"
#include <unordered_map>

int main() {
  IntervalTree test_amplicons = IntervalTree();

  test_amplicons.insert(Interval(1, 5));
  test_amplicons.insert(Interval(5, 30));

  std::vector<ITNode*> amp1;
  test_amplicons.find_read_amplicon(1, 2, amp1);

  std::vector<ITNode*> amp2;
  test_amplicons.find_read_amplicon(20, 25, amp2);

  site_aggregator sa;
  site_state ss1 = site_state();
  ss1.set_nucleotide("A", 30, 1);
  ss1.set_amplicon(amp1[0], false);
  site_state ss2 = site_state();
  ss2.set_nucleotide("T", 25, 1);
  ss2.set_amplicon(amp1[0], false);
  site_state ss3 = site_state();
  ss3.set_nucleotide("T", 45, 1);
  ss3.set_amplicon(amp1[0], false);

  site_state ss4 = site_state();
  ss4.set_nucleotide("-T", 22, 20);
  ss4.set_amplicon(amp1[0], false);
  site_state ss5 = site_state();
  ss5.set_nucleotide("C", 20, 20);
  ss5.set_amplicon(amp1[0], false);
  site_state ss6 = site_state();
  ss6.set_nucleotide("+G", 20, 20);
  ss6.set_amplicon(amp1[0], false);
  site_state ss7 = site_state();
  ss7.set_nucleotide("+G", 20, 20);
  ss7.set_amplicon(amp1[0], false);

  site_state ss8 = site_state();
  ss8.set_nucleotide("G", 15, 50);
  ss8.set_amplicon(nullptr, true);

  std::vector<site_state> site_states = {ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8};
  sa.aggregate(site_states);
  std::unordered_map<site_aggregator_key, site_aggregator_stats> aggregated_site_states = sa.get_data();
  for (const auto& [key, value] : aggregated_site_states) {
    std::cout << key.coordinate.position << " " << key.state;
    if(key.amplicon != nullptr)
      std::cout << " (" << key.amplicon->data->low << ", " << key.amplicon->data->high << ")" << std::endl;
    else
      std::cout << " (, )" << std::endl;
    std::cout << value.get_pos() << " " << value.get_depth() << " " << value.get_gapped_depth() << " " << value.get_qual() << std::endl << std::endl;
  }
  return 0;
}