#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

//prints mismatches in the same format as the checks below
bool check_base(const consensus_sequence &cs, uint32_t scenario, uint32_t genome, uint32_t pos, std::string expected){
  std::string actual = cs.get_base(pos);
  if(actual != expected){
    std::cerr << "scenario " << scenario << " genome" << genome << " position " << pos << " expected " << expected
              << " got " << actual << std::endl;
    return false;
  }
  return true;
}

int main() {
  int num_tests = 0;
  int success = 0;

  // FIRST SCENARIO (TESTs 1 and 2): 5 positions. Position 2 has three competing calls:
  //   -AAA  freq=0.90  (dominant deletion)  -> consensus genome 0
  //   -AA   freq=0.03  (minor, overlapped)  -> excluded (would be suppressed
  //                                            upstream by set_deletion_flags,
  //                                            so it never reaches either genome)
  //   A     freq=0.10  (ref persists)       -> consensus genome 1
  // All other positions (1,3,4,5) are background/invariant: a single "A" call,
  // half_normal_upper=true, assigned to both genomes.

  uint32_t max_position = 5;
  {
    variant background1;
    background1.position = 1;
    background1.nuc = "A";
    background1.half_normal_upper = true;
    background1.consensus_numbers = {0, 1}; //background, both genomes

    variant major_del;
    major_del.position = 2;
    major_del.nuc = "-AAA";
    major_del.consensus_numbers = {0}; //dominant deletion, freq 0.90, spans 2-4

    variant minor_ref;
    minor_ref.position = 2;
    minor_ref.nuc = "A";
    minor_ref.consensus_numbers = {1}; //reference persists, freq 0.10

    variant minor_del;
    minor_del.position = 2;
    minor_del.nuc = "-AA";
    minor_del.consensus_numbers = {1}; //deletion which is likely an error, freq 0.03, spans 2-3

    variant minor_ref3;
    minor_ref3.position = 3;
    minor_ref3.nuc = "A";
    minor_ref3.consensus_numbers = {1}; //genome 1 only, freq 0.10

    variant minor_ref4;
    minor_ref4.position = 4;
    minor_ref4.nuc = "A";
    minor_ref4.consensus_numbers = {1}; //genome 1 only, freq 0.10

    variant background5;
    background5.position = 5;
    background5.nuc = "A";
    background5.half_normal_upper = true;
    background5.consensus_numbers = {0, 1}; //background, both genomes

    std::vector<variant> scenario1 = {background1, major_del, minor_ref, minor_del, minor_ref3, minor_ref4, background5};

    std::vector<consensus_sequence> genomes;
    for(uint32_t i=0; i < 2; i++){
      genomes.emplace_back(max_position);
    }
    assign_variants_position(scenario1, genomes);
    for(uint32_t i=0; i < genomes.size(); i++){
      genomes[i].process_variant_assignments();
      genomes[i].get_consensus();
    }
    consensus_sequence &genome0 = genomes[0];
    consensus_sequence &genome1 = genomes[1];

    // TEST 1 - genome 0 (the dominant/majority genome) should show the deletion
    // at position 2 and "A" everywhere else.
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = (pos == 2 || pos == 3 || pos == 4) ? "-" : "A";
        std::string actual = genome0.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "scenario 1 genome0 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }

    // TEST 2 - genome 1 (the minor genome, no deletion) should show "A" at every
    // position, except for 2 and 3 which are 'N'
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = "A";
        if (pos == 2 || pos == 3) expected = "N";
        std::string actual = genome1.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "scenario 1 genome1 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }
  }

  // SECOND SCENARIO (TESTs 3 and 4): 5 positions. Position 2 has two competing alleles,
  // no deletion this time: A@0.90 -> genome 0, T@0.10 -> genome 1, each its own
  // standalone record. Position 3 has a background invariant "A"
  // (half_normal_upper, both genomes) plus a minor insertion +CC@0.10 assigned
  // only to genome 1.
  // All other positions are background "A"s (half_normal_upper,
  // both genomes).
  {
    variant background1;
    background1.position = 1;
    background1.nuc = "A";
    background1.half_normal_upper = true;
    background1.consensus_numbers = {0, 1}; //background, both genomes

    variant major_allele;
    major_allele.position = 2;
    major_allele.nuc = "A";
    major_allele.consensus_numbers = {0}; //freq 0.90

    variant minor_allele;
    minor_allele.position = 2;
    minor_allele.nuc = "T";
    minor_allele.consensus_numbers = {1}; //freq 0.10

    variant background3;
    background3.position = 3;
    background3.nuc = "A";
    background3.half_normal_upper = true;
    background3.consensus_numbers = {0, 1}; //background, both genomes

    variant minor_insertion;
    minor_insertion.position = 3;
    minor_insertion.nuc = "+CC";
    minor_insertion.consensus_numbers = {1}; //freq 0.10

    variant background4;
    background4.position = 4;
    background4.nuc = "A";
    background4.half_normal_upper = true;
    background4.consensus_numbers = {0, 1}; //background, both genomes

    variant background5;
    background5.position = 5;
    background5.nuc = "A";
    background5.half_normal_upper = true;
    background5.consensus_numbers = {0, 1}; //background, both genomes

    std::vector<variant> scenario2 = {background1, major_allele, minor_allele, background3, minor_insertion, background4, background5};

    std::vector<consensus_sequence> genomes;
    for(uint32_t i=0; i < 2; i++){
      genomes.emplace_back(max_position);
    }
    assign_variants_position(scenario2, genomes);
    for(uint32_t i=0; i < genomes.size(); i++){
      genomes[i].process_variant_assignments();
      genomes[i].get_consensus();
    }
    consensus_sequence &genome0 = genomes[0];
    consensus_sequence &genome1 = genomes[1];

    // TEST 3 - genome 0: plain "A" at every position (no insertion, no minor allele).
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = "A";
        std::string actual = genome0.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "scenario 2 genome0 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }

    // TEST 4 - genome 1: "T" at position 2 (the minor allele), "ACC" at
    // position 3 (background "A" with the insertion appended), "A" elsewhere.
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = "A";
        if (pos == 2) expected = "T";
        if (pos == 3) expected = "ACC";
        std::string actual = genome1.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "scenario 2 genome1 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }
  }

  // THIRD SCENARIO (TESTs 5 and 6): ambiguities. Starting baseline - 5 positions, all
  // background/invariant "A" (half_normal_upper=true, both genomes), same
  // style as the first two scenarios' background positions.
  {
    variant background1;
    background1.position = 1;
    background1.nuc = "A";
    background1.half_normal_upper = true;
    background1.consensus_numbers = {0, 1}; //background, both genomes

    variant major_allele;
    major_allele.position = 2;
    major_allele.nuc = "A";
    major_allele.consensus_numbers = {0}; //freq 0.90

    variant minor_c;
    minor_c.position = 2;
    minor_c.nuc = "C";
    minor_c.consensus_numbers = {1}; //freq 0.06

    variant minor_t;
    minor_t.position = 2;
    minor_t.nuc = "T";
    minor_t.consensus_numbers = {1}; //freq 0.04

    variant background3;
    background3.position = 3;
    background3.nuc = "A";
    background3.half_normal_upper = true;
    background3.consensus_numbers = {0, 1}; //background, both genomes

    variant vague_a;
    vague_a.position = 4;
    vague_a.nuc = "A";
    vague_a.consensus_numbers = {0, 1}; //vague assignment, freq 0.66

    variant vague_c;
    vague_c.position = 4;
    vague_c.nuc = "C";
    vague_c.consensus_numbers = {0, 1}; //vague assignment, freq 0.44

    variant background5;
    background5.position = 5;
    background5.nuc = "A";
    background5.half_normal_upper = true;
    background5.consensus_numbers = {0, 1}; //background, both genomes

    std::vector<variant> scenario3 = {background1, major_allele, minor_c, minor_t, background3, vague_a, vague_c, background5};

    std::vector<consensus_sequence> genomes;
    for(uint32_t i=0; i < 2; i++){
      genomes.emplace_back(max_position);
    }
    assign_variants_position(scenario3, genomes);
    for(uint32_t i=0; i < genomes.size(); i++){
      genomes[i].process_variant_assignments();
      genomes[i].get_consensus();
    }
    consensus_sequence &genome0 = genomes[0];
    consensus_sequence &genome1 = genomes[1];

    // TEST 5 - genome 0: "A" at every position except position 4, where the
    // vaguely-assigned "A"/"C" alleles should combine into IUPAC code "M".
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = (pos == 4) ? "M" : "A";
        std::string actual = genome0.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "scenario 3 genome0 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }

    // TEST 6 - genome 1: "A" everywhere except position 2, where the
    // competing "C"/"T" minor alleles should combine into the IUPAC
    // ambiguity code "Y", and position 4, where the vaguely-assigned
    // "A"/"C" alleles should combine into IUPAC code "M".
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = "A";
        if (pos == 2) expected = "Y";
        if (pos == 4) expected = "M";
        std::string actual = genome1.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "scenario 3 genome1 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }
  }

  // FOURTH SCENARIO (TESTs 7 and 8) - a genome listed in a variant's ambiguous_numbers is N at that position,
  // even when other variants confidently assign it a base there
  {
    variant v_confirmed;
    v_confirmed.position = 5;
    v_confirmed.nuc = "A";
    v_confirmed.gapped_freq = 0.40;
    v_confirmed.consensus_numbers = {3}; 

    variant v_ambiguous;
    v_ambiguous.position = 5;
    v_ambiguous.nuc = "T";
    v_ambiguous.gapped_freq = 0.50;
    v_ambiguous.consensus_numbers = {1, 3};
    v_ambiguous.ambiguous_numbers = {1, 3}; 

    variant v_ambiguous2;
    v_ambiguous2.position = 5;
    v_ambiguous2.nuc = "G";
    v_ambiguous2.gapped_freq = 0.10;
    v_ambiguous2.consensus_numbers = {1};

    variant v_unambiguous;
    v_unambiguous.position = 5;
    v_unambiguous.nuc = "C";
    v_unambiguous.gapped_freq = 0.10;
    v_unambiguous.consensus_numbers = {0};

    std::vector<variant> ambiguity_variants = {v_confirmed, v_ambiguous, v_ambiguous2, v_unambiguous};
    std::vector<consensus_sequence> ambiguity_consensus_seqs;
    for(uint32_t i=0; i < 4; i++){
      ambiguity_consensus_seqs.emplace_back(10);
    }
    assign_variants_position(ambiguity_variants, ambiguity_consensus_seqs);
    for(uint32_t i=0; i < ambiguity_consensus_seqs.size(); i++){
      ambiguity_consensus_seqs[i].process_variant_assignments();
      ambiguity_consensus_seqs[i].get_consensus();
    }

    //TEST 7 - check that base 5 is ambiguous in both genomes
    {
      bool pass = true;
      pass &= check_base(ambiguity_consensus_seqs[3], 4, 3, 5, "N");
      pass &= check_base(ambiguity_consensus_seqs[1], 4, 1, 5, "N");
      num_tests++;
      if (pass) success++;
    }

    //TEST 8 - check that genome 0 has the unambiguous base
    {
      bool pass = check_base(ambiguity_consensus_seqs[0], 4, 0, 5, "C");
      num_tests++;
      if (pass) success++;
    }
  }

  // FIFTH SCENARIO (TEST 9) - a deletion assigned to the upper half normal is written as a gap,
  // not as its deleted bases repeated once per span position
  {
    variant background1;
    background1.position = 1;
    background1.nuc = "G";
    background1.half_normal_upper = true;
    background1.consensus_numbers = {0};

    variant upper_del;
    upper_del.position = 2;
    upper_del.nuc = "-TAC"; //spans positions 2-4
    upper_del.half_normal_upper = true;
    upper_del.consensus_numbers = {0};

    variant background5;
    background5.position = 5;
    background5.nuc = "G";
    background5.half_normal_upper = true;
    background5.consensus_numbers = {0};

    variant background6;
    background6.position = 6;
    background6.nuc = "T";
    background6.half_normal_upper = true;
    background6.consensus_numbers = {0};

    variant background7;
    background7.position = 7;
    background7.nuc = "A";
    background7.half_normal_upper = true;
    background7.consensus_numbers = {0};

    variant background8;
    background8.position = 8;
    background8.nuc = "C";
    background8.half_normal_upper = true;
    background8.consensus_numbers = {0};
    std::vector<variant> upper_del_variants = {background1, upper_del, background5, background6,
                                               background7, background8};

    std::vector<consensus_sequence> upper_del_seqs;
    upper_del_seqs.emplace_back(8);
    assign_variants_position(upper_del_variants, upper_del_seqs);
    upper_del_seqs[0].process_variant_assignments();
    upper_del_seqs[0].get_consensus();

    //TEST 9 - check the deletion is written as a gap
    {
      std::string expected = "G---GTAC";
      bool pass = true;
      for(uint32_t pos = 1; pos <= expected.size(); pos++){
        pass &= check_base(upper_del_seqs[0], 5, 0, pos, std::string(1, expected[pos-1]));
      }
      num_tests++;
      if (pass) success++;
    }
  }

  // SIXTH SCENARIO (TESTs 10, 11 and 12) - a freq_outlier variant sits too far from the cluster it was assigned to for any
  // population to credibly carry it, so its position is N in every genome, including ones
  // it was never assigned to. A deletion Ns its whole span.
  {
    variant v_wide;
    v_wide.position = 1;
    v_wide.nuc = "A";
    v_wide.consensus_numbers = {0}; //assigned to genome 0 only
    v_wide.freq_outlier = true;

    variant v_ok;
    v_ok.position = 2;
    v_ok.nuc = "G";
    v_ok.consensus_numbers = {0}; //fits its cluster, must keep its base

    variant v_wide_del;
    v_wide_del.position = 3;
    v_wide_del.nuc = "-TA"; //spans positions 3 and 4
    v_wide_del.consensus_numbers = {0};
    v_wide_del.freq_outlier = true;

    variant v_ok2;
    v_ok2.position = 5;
    v_ok2.nuc = "C";
    v_ok2.consensus_numbers = {1};

    variant v_wide_ins;
    v_wide_ins.position = 5;
    v_wide_ins.nuc = "+CC";
    v_wide_ins.consensus_numbers = {1}; //in a different genomes because insertions can't follow deletions
    v_wide_ins.freq_outlier = true;

    std::vector<variant> wide_variants = {v_wide, v_ok, v_wide_del, v_ok2, v_wide_ins};
    std::vector<consensus_sequence> wide_seqs;
    for(uint32_t i=0; i < 2; i++){
      wide_seqs.emplace_back(10);
    }
    assign_variants_position(wide_variants, wide_seqs);
    for(uint32_t i=0; i < wide_seqs.size(); i++){
      wide_seqs[i].process_variant_assignments();
      wide_seqs[i].get_consensus();
    }

    //TEST 10 - the outlier substitution at 1 is N in every genome, and the variant at 2 keeps its base
    {
      bool pass = true;
      for(uint32_t i=0; i < wide_seqs.size(); i++){
        pass &= check_base(wide_seqs[i], 6, i, 1, "N");
      }
      pass &= check_base(wide_seqs[0], 6, 0, 2, "G");
      num_tests++;
      if (pass) success++;
    }

    //TEST 11 - the outlier deletion Ns its whole span, 3 and 4, in every genome
    {
      bool pass = true;
      for(uint32_t i=0; i < wide_seqs.size(); i++){
        pass &= check_base(wide_seqs[i], 6, i, 3, "N");
        pass &= check_base(wide_seqs[i], 6, i, 4, "N");
      }
      num_tests++;
      if (pass) success++;
    }

    //TEST 12 - the outlier insertion is dropped, and genome 1 keeps the base it follows at 5
    {
      bool pass = check_base(wide_seqs[1], 6, 1, 5, "C");
      num_tests++;
      if (pass) success++;
    }
  }

  // SEVENTH SCENARIO (TESTs 13, 14, 15 and 16) - a low quality or low depth call is an ambiguity in every genome it was
  // assigned to, no matter what: even in the upper half normal and across a deletion's span.
  // A low depth insertion is simply left out
  {
    variant v_low_depth;
    v_low_depth.position = 1;
    v_low_depth.nuc = "A";
    v_low_depth.half_normal_upper = true;
    v_low_depth.depth_flag = true;
    v_low_depth.consensus_numbers = {0, 1};

    variant v_low_qual;
    v_low_qual.position = 2;
    v_low_qual.nuc = "C";
    v_low_qual.half_normal_upper = true;
    v_low_qual.qual_flag = true;
    v_low_qual.consensus_numbers = {0, 1};

    variant v_ok;
    v_ok.position = 3;
    v_ok.nuc = "G";
    v_ok.half_normal_upper = true;
    v_ok.consensus_numbers = {0, 1}; //unflagged, must keep its base

    variant v_ref4;
    v_ref4.position = 4;
    v_ref4.nuc = "T";
    v_ref4.consensus_numbers = {0, 1};

    variant v_ref5;
    v_ref5.position = 5;
    v_ref5.nuc = "A";
    v_ref5.consensus_numbers = {0, 1};

    variant v_low_depth_del;
    v_low_depth_del.position = 4;
    v_low_depth_del.nuc = "-TA"; //spans positions 4 and 5
    v_low_depth_del.depth_flag = true;
    v_low_depth_del.consensus_numbers = {0};

    variant v_anchor;
    v_anchor.position = 6;
    v_anchor.nuc = "C";
    v_anchor.consensus_numbers = {0, 1};

    variant v_low_depth_ins;
    v_low_depth_ins.position = 6;
    v_low_depth_ins.nuc = "+GG";
    v_low_depth_ins.depth_flag = true;
    v_low_depth_ins.consensus_numbers = {1};

    variant v_lower;
    v_lower.position = 7;
    v_lower.nuc = "T";
    v_lower.half_normal_lower = true;
    v_lower.consensus_numbers = {0}; //lower half normal noise, must be ignored

    variant v_ref7;
    v_ref7.position = 7;
    v_ref7.nuc = "C";
    v_ref7.consensus_numbers = {0}; //normal call at 7 in genome 0

    std::vector<variant> flagged_variants = {v_low_depth, v_low_qual, v_ok, v_ref4, v_ref5,
                                             v_low_depth_del, v_anchor, v_low_depth_ins, v_lower, v_ref7};
    std::vector<consensus_sequence> flagged_seqs;
    for(uint32_t i=0; i < 2; i++){
      flagged_seqs.emplace_back(10);
    }
    assign_variants_position(flagged_variants, flagged_seqs);
    for(uint32_t i=0; i < flagged_seqs.size(); i++){
      flagged_seqs[i].process_variant_assignments();
      flagged_seqs[i].get_consensus();
    }

    //TEST 13 - low depth and low quality calls in the upper half normal are N in both genomes,
    //and the unflagged call at 3 keeps its base
    {
      bool pass = true;
      for(uint32_t i=0; i < flagged_seqs.size(); i++){
        pass &= check_base(flagged_seqs[i], 7, i, 1, "N");
        pass &= check_base(flagged_seqs[i], 7, i, 2, "N");
        pass &= check_base(flagged_seqs[i], 7, i, 3, "G");
      }
      num_tests++;
      if (pass) success++;
    }

    //TEST 14 - the low depth deletion Ns its whole span in genome 0, and genome 1 keeps its bases
    {
      bool pass = true;
      pass &= check_base(flagged_seqs[0], 7, 0, 4, "N");
      pass &= check_base(flagged_seqs[0], 7, 0, 5, "N");
      pass &= check_base(flagged_seqs[1], 7, 1, 4, "T");
      pass &= check_base(flagged_seqs[1], 7, 1, 5, "A");
      num_tests++;
      if (pass) success++;
    }

    //TEST 15 - the low depth insertion is left out, so both genomes keep the base it follows
    {
      bool pass = true;
      pass &= check_base(flagged_seqs[1], 7, 1, 6, "C");
      pass &= check_base(flagged_seqs[0], 7, 0, 6, "C");
      num_tests++;
      if (pass) success++;
    }

    //TEST 16 - the lower half normal call is ignored, so genome 0 keeps the normal call at 7
    {
      bool pass = check_base(flagged_seqs[0], 7, 0, 7, "C");
      num_tests++;
      if (pass) success++;
    }
  }

  // EIGHTH SCENARIO (TESTs 17 and 18) - a variant whose frequency fluctuates across amplicons (position_masked) is N
  // in the genome it was assigned to, while a different variant at the same position in another
  // genome keeps its base
  {
    variant v_masked;
    v_masked.position = 3;
    v_masked.nuc = "A";
    v_masked.position_masked = true;
    v_masked.consensus_numbers = {0};

    variant v_kept;
    v_kept.position = 3;
    v_kept.nuc = "G";
    v_kept.consensus_numbers = {1};

    std::vector<variant> masked_variants = {v_masked, v_kept};
    std::vector<consensus_sequence> masked_seqs;
    for(uint32_t i=0; i < 2; i++){
      masked_seqs.emplace_back(10);
    }
    assign_variants_position(masked_variants, masked_seqs);
    for(uint32_t i=0; i < masked_seqs.size(); i++){
      masked_seqs[i].process_variant_assignments();
      masked_seqs[i].get_consensus();
    }

    //TEST 17 - the position masked variant is N in genome 0
    {
      bool pass = check_base(masked_seqs[0], 8, 0, 3, "N");
      num_tests++;
      if (pass) success++;
    }

    //TEST 18 - the unmasked variant at the same position keeps its base in genome 1
    {
      bool pass = check_base(masked_seqs[1], 8, 1, 3, "G");
      num_tests++;
      if (pass) success++;
    }
  }

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
