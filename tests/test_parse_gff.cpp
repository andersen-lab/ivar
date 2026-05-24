#include "../src/parse_gff.h"

int main() {
  int success = 0;
  gff3 g("../../data/test_discontinuous_cds.gff");

  const std::vector<cds_group> &groups = g.get_cds_groups();

  success += groups.size() == 3;
  if(success != 1) {
    std::cout << "Error in number of groups";
    return -1;
  }


  // Forward discontinuous
  const cds_group *fwd = g.find_cds_group_by_id("id-fwd-discontinuous");
  success += fwd != nullptr;
  if(success != 2) {
    std::cout << "Error in finding fwd group";
    return -1;
  }

  success += fwd->segments().size() == 2;
  if(success != 3) {
    std::cout << "Error in number of segments in fwd group";
    return -1;
  }

  success += (fwd->segments()[0].get_start() == 10) && (fwd->segments()[0].get_end() == 29);
  if(success != 4) {
    std::cout << "Error in coordinates of start and end of first segment in fwd group";
    return -1;
  }

  success += (fwd->segments()[1].get_start() == 100) && (fwd->segments()[1].get_end() == 149);
  if(success != 5) {
    std::cout << "Error in coordinates of start and end of second segment in fwd group";
    return -1;
  }

  success += (fwd->get_strand() == '+') && (fwd->get_phase() == 0) && (fwd->length() == 70);
  if(success != 6) {
    std::cout << "Error in strand, phase of length of group";
    return -1;
  }

  int64_t contains_nuc_pos[4] = {15, 50, 120, 150};
  bool expected_contains[4] = {true, false, true, false};

  for(int i = 0; i < 4; i++) {
    success += (fwd->contains(contains_nuc_pos[i]) == expected_contains[i]);
    if(success != 7 + i) {
      std::cout << "Error in contains for position " << contains_nuc_pos[i] << std::endl;
    }
  }

  int64_t nuc_pos_to_cds[5] = {10, 29, 100, 149, 50};
  int64_t expected_cds_pos[5] = {0, 19, 20, 69, -1};

  for(int i =0; i < 5; i++) {
    success += (fwd->genomic_to_cds_pos(nuc_pos_to_cds[i]) == expected_cds_pos[i]);
    if(success != 11 + i) {
      std::cout << "Error in genomic to cds position for position " << nuc_pos_to_cds[i] << std::endl;
    }
  }

  int64_t cds_pos_to_nuc[5] = {0, 19, 20, 69, 70};
  int64_t expected_nuc_pos[5] = {10, 29, 100, 149, -1};

  for (int i = 0; i < 5; ++i) {
    success += (fwd->cds_to_genomic_pos(cds_pos_to_nuc[i]) == expected_nuc_pos[i]);
    if(success != 16 + i) {
      std::cout << "Error in cds to genomic position for position " << cds_pos_to_nuc[i] << std::endl;
    }
  }

  // Reverse discontinuous segment
  const cds_group *rev = g.find_cds_group_by_id("id-rev-discontinuous");
  success += rev != nullptr;
  if(success != 21) {
    std::cout << "(rev) Error in finding rev group";
    return -1;
  }

  success += rev->segments().size() == 2;
  if(success != 22) {
    std::cout << "(rev) Error in number of segments in rev group";
    return -1;
  }

  success += (rev->segments()[0].get_start() == 300) && (rev->segments()[0].get_end() == 319);
  if(success != 23) {
    std::cout << "(rev) Error in coordinates of start and end of first segment in rev group";
    return -1;
  }

  success += (rev->segments()[1].get_start() == 200) && (rev->segments()[1].get_end() == 220);
  if(success != 24) {
    std::cout << "(rev) Error in coordinates of start and end of second segment in rev group";
    return -1;
  }

  success += (rev->get_strand() == '-') && (rev->get_phase() == 0) && (rev->length() == 41);
  if(success != 25) {
    std::cout << "(rev) Error in strand, phase of length of rev group";
    return -1;
  }

  int64_t nuc_pos_to_cds_rev[5] = {319, 300, 220, 200, 250};
  int64_t expected_cds_pos_rev[5] = {0, 19, 20, 40, -1};

  for(int i =0; i < 5; i++) {
    success += (rev->genomic_to_cds_pos(nuc_pos_to_cds_rev[i]) == expected_cds_pos_rev[i]);
    if(success != 26 + i) {
      std::cout << "(rev) Error in genomic to cds position for position " << nuc_pos_to_cds[i] << std::endl;
    }
  }

  int64_t cds_pos_to_nuc_rev[5] = {0, 19, 20, 40, 70};
  int64_t expected_nuc_pos_rev[5] = {319, 300, 220, 200, -1};

  for (int i = 0; i < 5; ++i) {
    success += (rev->cds_to_genomic_pos(cds_pos_to_nuc_rev[i]) == expected_nuc_pos_rev[i]);
    if(success != 31 + i) {
      std::cout << "(rev) Error in cds to genomic position for position " << cds_pos_to_nuc_rev[i] << std::endl;
    }
  }

  // query_cds_groups
  success += (g.query_cds_groups(15).size() == 1);
  success += (g.query_cds_groups(60).empty());

  return success == 37 ? 0 : 1;
}
