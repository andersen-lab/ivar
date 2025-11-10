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
  ss.set_nucleotide_gap(min_qual, 24);
  expected_read_site_states[2].push_back(ss);

  // Read name: test_20_381_1:0:0_6:0:0_1_150
  ss.set_nucleotide("-ACCA", 20, 353);
  expected_read_site_states[3].push_back(ss);
  ss.set_nucleotide_gap(min_qual, 354);
  expected_read_site_states[3].push_back(ss);
  ss.set_nucleotide_gap(min_qual, 355);
  expected_read_site_states[3].push_back(ss);
  ss.set_nucleotide_gap(min_qual, 356);
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
        if(expected_read_site_states[test_idx][j] == read_site_states[i]) {
          pass += 1;
        }
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

  // Test merge_reads()
  std::vector<site_state> read_site_states_one;
  std::vector<site_state> read_site_states_two;
  std::vector<site_state> merged_site_states;

  // Read 1: A  C (+TC) G (-CT) T  G +GT C
  // Qual:  20 20   10 30 20   30 40 20  30
  std::vector<std::string> read1_nuc = {"A", "C", "+TC", "G", "-CT", "T", "G", "+GT", "C"};
  std::vector<uint8_t> read1_qual = {20, 10, 20, 30, 20, 30, 40, 20, 20};
  for(int i = 0; i < read1_nuc.size(); i++) {
    site_state ss1;
    ss1.set_nucleotide(read1_nuc[i], read1_qual[i], i+1);
    read_site_states_one.push_back(ss1);
  }

  // Read 2:    C (+GT)  G (-CT) G  G +TG C
  // Qual:     10   10   30 20   30 20 20  20
  // Read starts at position 2
  std::vector<std::string> read2_nuc = {"C", "+GT", "G", "-CT", "G", "G", "+GT", "C"};
  std::vector<uint8_t> read2_qual = {10, 20, 30, 20, 30, 20, 20, 20};
  for(int i = 0; i < read2_nuc.size(); i++) {
    site_state ss2;
    ss2.set_nucleotide(read2_nuc[i], read2_qual[i], i+2);
    read_site_states_two.push_back(ss2);
  }

  vc.merge_reads(read_site_states_one, read_site_states_two, merged_site_states);

  std::vector<uint32_t> expected_merged_positions = {1, 2, 4, 5, 7, 8, 9};
  std::vector<std::string> expected_merged_nuc = {"A", "C", "G", "-CT", "G", "+GT", "C"};
  std::vector<uint8_t> expected_merged_qual = {20, 10, 30, 20, 30, 20, 20};


  if(expected_merged_positions.size() != merged_site_states.size()) {
    std::cerr << "Merged site states size incorrect." << std::endl;
    return 1;
  }
  for(int i = 0; i < merged_site_states.size(); i++){
    if(merged_site_states[i].coordinate.position != expected_merged_positions[i] ||
        merged_site_states[i].state != expected_merged_nuc[i] ||
        merged_site_states[i].quality != expected_merged_qual[i]) {
      std::cerr << "Merged site states did not match at index " << i << std::endl;
      return 1;
    }
  }

  sam_hdr_destroy(header);
  sam_close(bam);

  return 0;
}