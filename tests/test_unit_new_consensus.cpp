#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include "htslib/sam.h"
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/ref_seq.h"
#include "../src/parse_gff.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"
#include "../src/interval_tree.h"

int main() {
  int num_tests = 6;
  int success = 0;
  std::vector<variant> variants;
  std::vector<double> means = {0.90, 0.10};


  /* TEST 1 - set_deletion_flags: minor deletion starting at the same position as a dominant
   * deletion is marked overlapped_deletion=true; dominant is not.
   */
  {
    variant del_dom, del_min;
    del_dom.nuc = "-GGA"; del_dom.position = 5; del_dom.gapped_freq = 0.90;
    del_min.nuc = "-G";   del_min.position = 5; del_min.gapped_freq = 0.05;
    variants.clear();
    variants.push_back(del_dom);
    variants.push_back(del_min);
    set_deletion_flags(variants, 0.03, 0.03);
    if (!variants[0].overlapped_deletion && variants[1].overlapped_deletion) success++;
    else std::cerr << "test 1 failed" << std::endl;
  }

  /* TEST 2 - set_deletion_flags: minor deletion with a different start but a range that
   * overlaps the dominant deletion is still marked overlapped_deletion=true.
   */
  {
    variant del_dom, del_min;
    del_dom.nuc = "-GGA"; del_dom.position = 5; del_dom.gapped_freq = 0.90;
    del_min.nuc = "-AG";  del_min.position = 6; del_min.gapped_freq = 0.05;
    variants.clear();
    variants.push_back(del_dom);
    variants.push_back(del_min);
    set_deletion_flags(variants, 0.03, 0.03);
    if (!variants[0].overlapped_deletion && variants[1].overlapped_deletion) success++;
    else std::cerr << "test 2 failed" << std::endl;
  }

  /* TEST 3 - set_deletion_flags: two deletions whose ranges do not overlap are both kept
   * (neither gets overlapped_deletion=true).
   */
  {
    variant del_a, del_b;
    del_a.nuc = "-G"; del_a.position = 5; del_a.gapped_freq = 0.90;
    del_b.nuc = "-G"; del_b.position = 7; del_b.gapped_freq = 0.05;
    variants.clear();
    variants.push_back(del_a);
    variants.push_back(del_b);
    set_deletion_flags(variants, 0.03, 0.03);
    if (!variants[0].overlapped_deletion && !variants[1].overlapped_deletion) success++;
    else std::cerr << "test 3 failed" << std::endl;
  }

  /* TEST 4 - flag_position_conflicts: two variants at the same position assigned to the
   * same cluster both receive position_conflict=true.
   */
  {
    variant va, vb;
    va.position = 10; va.nuc = "A"; va.cluster_assigned = 0;
    vb.position = 10; vb.nuc = "T"; vb.cluster_assigned = 0;
    variants.clear();
    variants.push_back(va);
    variants.push_back(vb);
    flag_position_conflicts(variants);
    if (variants[0].position_conflict && variants[1].position_conflict) success++;
    else std::cerr << "test 4 failed" << std::endl;
  }

  /* TEST 5 - flag_position_conflicts: two variants at the same position assigned to
   * different clusters do not get flagged.
   */
  {
    variant va, vb;
    va.position = 10; va.nuc = "A"; va.cluster_assigned = 0;
    vb.position = 10; vb.nuc = "T"; vb.cluster_assigned = 1;
    variants.clear();
    variants.push_back(va);
    variants.push_back(vb);
    flag_position_conflicts(variants);
    if (!variants[0].position_conflict && !variants[1].position_conflict) success++;
    else std::cerr << "test 5 failed" << std::endl;
  }

  /* TEST 6 - flag_position_conflicts: a deletion spans positions 5-7 (nuc="-GGA"). Its start
   * position (5) does not conflict with anything, but its last covered position (7) collides
   * with another variant assigned to the same cluster. The deletion, which never had its own
   * start flagged, must still pick up position_conflict via its span, and so must the variant
   * it collides with at position 7. The variant at position 5 (different cluster) must not.
   */
  {
    variant del, v5, v7;
    del.position = 5; del.nuc = "-GGA"; del.cluster_assigned = 0;
    v5.position = 5;  v5.nuc = "C"; v5.cluster_assigned = 1;
    v7.position = 7;  v7.nuc = "T"; v7.cluster_assigned = 0;
    variants.clear();
    variants.push_back(del);
    variants.push_back(v5);
    variants.push_back(v7);
    flag_position_conflicts(variants);
    if (variants[0].position_conflict && !variants[1].position_conflict && variants[2].position_conflict) success++;
    else std::cerr << "test 6 failed" << std::endl;
  }
  std::cerr << "success " << success << " tests " << num_tests << std::endl;
  return (num_tests == success) ? 0 : -1;
}
