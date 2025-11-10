#include "../src/variant_caller.h"

void print_cigar(uint32_t *cigar, int nlength) {
  for (int i = 0; i < nlength; ++i) {
    std::cerr << ((cigar[i]) & BAM_CIGAR_MASK);
    std::cerr << "-" << ((cigar[i]) >> BAM_CIGAR_SHIFT) << " ";
  }
  std::cerr << std::endl;
}

int main(){
  uint8_t min_qual = 20;
  std::string ref_path="../../data/db/test_ref.fa";
  variant_caller vc(min_qual, ref_path, "");
  samFile *bam = sam_open("../../data/test.unmapped.sorted.bam", "r");
  sam_hdr_t *header = sam_hdr_read(bam);
  bam1_t *aln = bam_init1();
  int ctr = 0, test_idx = 0;

  // TODO: Initialize expected_read_site_states better
  std::vector<std::vector<site_state>> expected_read_site_states(4);
  site_state ss;
  // Read name: test_outside_primer
  ss.set_nucleotide("-C", 20, 9);
  expected_read_site_states[0].push_back(ss);
  ss.set_nucleotide("G", 52, 42);
  expected_read_site_states[0].push_back(ss);
  // Read name: test_20_381_1:0:0_6:0:0_1_150
  ss.set_nucleotide("+G", 63, 30);
  expected_read_site_states[1].push_back(ss);
  ss.set_nucleotide("T", 49, 42);
  expected_read_site_states[1].push_back(ss);
  // Read name: test_20_381_1_softclip:0:0_6:0:0_1_150
  ss.set_nucleotide("C", 20, 20);
  expected_read_site_states[2].push_back(ss);
  ss.set_nucleotide("+A", 63, 25);
  expected_read_site_states[2].push_back(ss);
  ss.set_nucleotide("-CT", 20, 23);
  expected_read_site_states[2].push_back(ss);
  // Read name: test_20_381_1:0:0_6:0:0_1_150
  ss.set_nucleotide("-ACCA", 20, 353);
  expected_read_site_states[3].push_back(ss);
  ss.set_nucleotide("T", 27, 365);
  expected_read_site_states[3].push_back(ss);

  IntervalTree test_amplicons = IntervalTree();

  test_amplicons.insert(Interval(0, 199));
  test_amplicons.insert(Interval(10, 30));

  vc.set_amplicons(test_amplicons);

  // Note that amplicon assignment is based on interval_tree and find_read_amplicon() is tested there
  std::vector<bool> expected_ambiguous = {false, false, false, true};

  // Test parse_read()
  while (sam_read1(bam, header, aln) >= 0) {
    if (ctr != 0 && ctr != 1 && ctr != 2 && ctr != 7) {
      ctr ++;
      continue;
    }

    std::vector<site_state> read_site_states;
    vc.parse_read(aln, "test", read_site_states);
    std::vector<ITNode*> nodes;
    vc.assign_amplicon_to_read(aln->core.pos, bam_endpos(aln), read_site_states);


    int pass = 0;
    for(int i = 0; i < read_site_states.size(); i++){
      for (int j = 0; j < expected_read_site_states[test_idx].size(); j++){
        if(expected_read_site_states[test_idx][j] == read_site_states[i])
          pass += 1;
      }
    }

    if(pass != expected_read_site_states[test_idx].size()) {
      std::cerr << "Read parsing did not match for " << bam_get_qname(aln) << std::endl;
      return 1;
    }

    // Check amplicon assignment
    if(expected_ambiguous[test_idx] != read_site_states[0].is_ambiguous) {
      std::cerr << "Amplicon ambiguous flag did not match for " << bam_get_qname(aln) << std::endl;
      return 1;
    }

    ctr++;
    test_idx++;
  }

  // Test add_variants()
//  positions.assign({10, 11, 14, 17});
//  bases.assign({"A", "C", "-GT", "+A"});
//  qualities.assign({30, 40, min_qual, 60});
//
//  vc.add_variants(positions, bases, qualities);
//
//  positions.assign({11, 11, 14, 17});
//  bases.assign({"+AT", "T", "C", "C"});
//  qualities.assign({30, 40, min_qual, 60});
//
//  vc.add_variants(positions, bases, qualities);
//
//  std::vector<genomic_position> global_positions = vc.get_global_positions();
//
//  // TODO: Include expected depth, gapped_path, and alleles in a vector for testing
//  if(global_positions[10].depth != 1 ||
//      global_positions[10].gapped_depth != 1 ||
//      global_positions[10].alleles[0].nuc != "A")
//    return 1;
//
//  if(global_positions[14].depth != 1 ||
//      global_positions[14].gapped_depth != 2 ||
//      global_positions[14].alleles[0].nuc != "-GT")
//    return 1;
//
//  if(global_positions[11].depth != 2 ||
//      global_positions[11].gapped_depth != 2 ||
//      global_positions[11].alleles[1].nuc != "+AT")
//    return 1;

  return 0;
}