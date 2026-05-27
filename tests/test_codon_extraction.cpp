#include <iostream>
#include <string>
#include <vector>

#include "../src/site_state.h"
#include "../src/variant_caller.h"

static site_state make_nuc(char base, uint32_t pos, uint8_t qual = 30) {
  return site_state(base, qual, pos, NUCLEOTIDE);
}

int main() {
  variant_caller vc(20, "../../data/db/test_ref.fa", "../../data/test_discontinuous_cds.gff");
  vc.initialize_region("test");

  int success = 0;

  // Test 1: Ref: AGC (S). Alt ACT (T)
  std::vector<site_state> nuc_states;
  nuc_states.push_back(make_nuc('A', 10));
  nuc_states.push_back(make_nuc('C', 11));
  nuc_states.push_back(make_nuc('T', 12));
  std::vector<std::vector<site_state> > codon_by_grp, aa_by_grp;
  vc.extract_codon_and_aa_states(nuc_states, codon_by_grp, aa_by_grp);

  // 1 fwd (2 CDS), 1 rev (2 CDS), 1 fwd (1 CDS)
  if (codon_by_grp.size() == 3)
    success += 1;
  else
    std::cout << "Test 1: expected 3 groups, got " << codon_by_grp.size() << std::endl;

  if (codon_by_grp[0].size() == 1)
    success += 1;
  else
    std::cout << "Test 1: fwd group should have 1 codon, got " << codon_by_grp[0].size() << std::endl;

  if (codon_by_grp[0][0].state == "ACT")
    success += 1;
  else
    std::cout << "Test 1: codon state should be ACT, got '" << codon_by_grp[0][0].state << "'" << std::endl;

  if (aa_by_grp[0][0].state == "T")
    success += 1;
  else
    std::cout << "Test 1: AA should be T, got '" << aa_by_grp[0][0].state << "'" << std::endl;

  if (codon_by_grp[1].empty())
    success += 1;
  else
    std::cout << "Test 1: rev group should be empty, got " << codon_by_grp[1].size() << std::endl;

  if (codon_by_grp[2].empty())
    success += 1;
  else
    std::cout << "Test 1: single-segment group should be empty, got " << codon_by_grp[2].size() << std::endl;

  // Test 2: Only pos 10, 11 supplied — no codon emitted.
  nuc_states.clear();
  nuc_states.push_back(make_nuc('T', 10));
  nuc_states.push_back(make_nuc('G', 11));
  codon_by_grp.clear();
  aa_by_grp.clear();
  vc.extract_codon_and_aa_states(nuc_states, codon_by_grp, aa_by_grp);

  if (codon_by_grp[0].empty())
    success += 1;
  else
    std::cout << "Test 2: should be empty on partial cover, got " << codon_by_grp[0].size() << std::endl;

  // Test 3: Reverse-strand codon. Ref: GAV (V) 319,318,317. Alt: TAA
  nuc_states.clear();
  nuc_states.push_back(make_nuc('A', 317));
  nuc_states.push_back(make_nuc('A', 318));
  nuc_states.push_back(make_nuc('T', 319));
  codon_by_grp.clear();
  aa_by_grp.clear();
  vc.extract_codon_and_aa_states(nuc_states, codon_by_grp, aa_by_grp);

  if (codon_by_grp[1].size() == 1)
    success += 1;
  else
    std::cout << "Test 3: rev group should have 1 codon, got " << codon_by_grp[1].size() << std::endl;

  if (codon_by_grp[1][0].state == "ATT")
    success += 1;
  else
    std::cout << "Test 3: rev codon should be ATT, got '" << codon_by_grp[1][0].state << "'" << std::endl;

  if (aa_by_grp[1][0].state == "I")
    success += 1;
  else
    std::cout << "Test 3: rev AA should be I, got '" << aa_by_grp[1][0].state << "'" << std::endl;

  // Test 4: Codon across two CDS. Ref codon TCT (A). Alt codon TCA (S).
  nuc_states.clear();
  nuc_states.push_back(make_nuc('T', 28));
  nuc_states.push_back(make_nuc('C', 29));
  nuc_states.push_back(make_nuc('A', 100));
  codon_by_grp.clear();
  aa_by_grp.clear();
  vc.extract_codon_and_aa_states(nuc_states, codon_by_grp, aa_by_grp);

  bool found = false;
  // Sixth codon in ref
  if(codon_by_grp[0][0].coordinate.position == 6)
    success += 1;
  else
    std::cout << "Test 4: codon should be at position 6, got " << codon_by_grp[0][0].coordinate.position << std::endl;

  if (codon_by_grp[0][0].state == "TCA")
    success += 1;
  else
    std::cout << "Test 4: codon should be TCA, got '" << codon_by_grp[0][0].state << "'" << std::endl;

  if (aa_by_grp[0][0].state == "S")
    success += 1;
  else
    std::cout << "Test 4: AA should be S, got '" << aa_by_grp[0][0].state << "'" << std::endl;

  return success == 13 ? 0 : 1;
}
