#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

int main() {
  int num_tests = 0;
  int success = 0;

  // Unit test for consensus_sequence::get_consensus, built by hand (bypassing
  // gmm_model/assign_variants_position/process_variant_assignments entirely) so
  // variant_records reflect exactly the scenario we want to check.
  //
  // 5 positions. Position 2 has three competing calls:
  //   -AAA  freq=0.90  (dominant deletion)  -> consensus genome 0
  //   -AA   freq=0.03  (minor, overlapped)  -> excluded (would be suppressed
  //                                            upstream by set_deletion_flags,
  //                                            so it never reaches either genome)
  //   A     freq=0.10  (ref persists)       -> consensus genome 1
  // All other positions (1,3,4,5) are background/invariant: a single "A" call,
  // half_normal_upper=true, assigned to both genomes.

  uint32_t max_position = 5;
  consensus_sequence genome0(max_position);
  consensus_sequence genome1(max_position);

  for (uint32_t pos = 1; pos <= max_position; pos++) {
    if (pos == 2) continue;
    variant v{};
    v.position = pos;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};

    genome0.add_variant(pos, v);
    genome1.add_variant(pos, v);
  }

  {
    variant v{};
    v.position = 2;
    v.nuc = "-AAA";
    v.gapped_freq = 0.90;
    v.assigned_deletion = true; // resolved deletion for this position
    v.consensus_numbers = {0};
    genome0.add_variant(2, v);
  }

  // -AA (freq 0.03) intentionally not added anywhere: represents the minor,
  // overlapped deletion that set_deletion_flags would suppress upstream.

  {
    variant v{};
    v.position = 2;
    v.nuc = "A";
    v.gapped_freq = 0.10;
    v.consensus_numbers = {1};
    genome1.add_variant(2, v);
  }

  genome0.get_consensus(0);
  genome1.get_consensus(1);

  // TEST 1 - genome 0 (the dominant/majority genome) should show the deletion
  // at position 2 and "A" everywhere else.
  {
    bool pass = true;
    for (uint32_t pos = 1; pos <= max_position; pos++) {
      std::string expected = (pos == 2) ? "-" : "A";
      std::string actual = genome0.get_base(pos);
      if (actual != expected) {
        pass = false;
        std::cerr << "genome0 position " << pos << " expected " << expected
                  << " got " << actual << std::endl;
      }
    }
    num_tests++;
    if (pass) success++;
  }

  // TEST 2 - genome 1 (the minor genome, no deletion) should show "A" at every
  // position, including position 2 where the reference persists.
  {
    bool pass = true;
    for (uint32_t pos = 1; pos <= max_position; pos++) {
      std::string expected = "A";
      std::string actual = genome1.get_base(pos);
      if (actual != expected) {
        pass = false;
        std::cerr << "genome1 position " << pos << " expected " << expected
                  << " got " << actual << std::endl;
      }
    }
    num_tests++;
    if (pass) success++;
  }

  // Second scenario: 5 positions again. Position 2 has two competing alleles,
  // no deletion this time: A@0.90 -> genome 0, T@0.10 -> genome 1, each its own
  // standalone record. Position 3 has a background invariant "A"
  // (half_normal_upper, both genomes) plus a minor insertion +CC@0.10 assigned
  // only to genome 1. All other positions are background "A"s (half_normal_upper,
  // both genomes).
  {
    consensus_sequence genome0(max_position);
    consensus_sequence genome1(max_position);

    for (uint32_t pos = 1; pos <= max_position; pos++) {
      if (pos == 2) continue;
      variant v{};
      v.position = pos;
      v.nuc = "A";
      v.gapped_freq = 1.0;
      v.half_normal_upper = true;
      v.position_half_normal_upper = true;
      v.consensus_numbers = {0, 1};

      genome0.add_variant(pos, v);
      genome1.add_variant(pos, v);
    }

    {
      variant v{};
      v.position = 2;
      v.nuc = "A";
      v.gapped_freq = 0.90;
      v.consensus_numbers = {0};
      genome0.add_variant(2, v);
    }

    {
      variant v{};
      v.position = 2;
      v.nuc = "T";
      v.gapped_freq = 0.10;
      v.consensus_numbers = {1};
      genome1.add_variant(2, v);
    }

    {
      variant v{};
      v.position = 3;
      v.nuc = "+CC";
      v.gapped_freq = 0.10;
      v.position_half_normal_upper = true; // position 3 also has the background half_normal_upper "A"
      v.consensus_numbers = {1};
      genome1.add_variant(3, v);
    }

    genome0.get_consensus(2);
    genome1.get_consensus(3);

    // TEST 3 - genome 0: plain "A" at every position (no insertion, no minor allele).
    {
      bool pass = true;
      for (uint32_t pos = 1; pos <= max_position; pos++) {
        std::string expected = "A";
        std::string actual = genome0.get_base(pos);
        if (actual != expected) {
          pass = false;
          std::cerr << "genome0 position " << pos << " expected " << expected
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
          std::cerr << "genome1 position " << pos << " expected " << expected
                    << " got " << actual << std::endl;
        }
      }
      num_tests++;
      if (pass) success++;
    }
  }

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
