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

  // FIRST SCENARIO: 5 positions. Position 2 has three competing calls:
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

  // genome 0, position 1: background "A"
  {
    variant v{};
    v.position = 1;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome0.add_variant(1, v);
  }
  // genome 1, position 1: background "A"
  {
    variant v{};
    v.position = 1;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome1.add_variant(1, v);
  }

  // genome 0, position 2: dominant deletion -AAA (freq 0.90), spans 2-4
  {
    variant v{};
    v.position = 2;
    v.nuc = "-AAA";
    v.gapped_freq = 0.90;
    v.assigned_deletion = true;
    v.consensus_numbers = {0};
    genome0.add_variant(2, v);
  }

  // genome 1, position 2 minor variant
  {
    variant v{};
    v.position = 2;
    v.nuc = "A";
    v.gapped_freq = 0.10;
    v.assigned_deletion = true;
    v.consensus_numbers = {1};
    genome1.add_variant(2, v);
  }
  // genome 1, position 2 deletion which is likely an error
  {
    variant v{};
    v.position = 2;
    v.nuc = "-AA";
    v.gapped_freq = 0.03;
    v.assigned_deletion = true;
    v.consensus_numbers = {1};
    genome1.add_variant(2, v);
  }
  // genome 0, position 3 doesn't exist
  // genome 1, position 3
  {
    variant v{};
    v.position = 3;
    v.nuc = "A";
    v.gapped_freq = 0.10;
    v.assigned_deletion = true;
    v.consensus_numbers = {1};
    genome1.add_variant(3, v);
  }
  // genome 0, position 4 doesn't exist
  // genome 1, position 4 deletion which is likely an error
  {
    variant v{};
    v.position = 4;
    v.nuc = "A";
    v.gapped_freq = 0.10;
    v.assigned_deletion = false;
    v.consensus_numbers = {1};
    genome1.add_variant(4, v);
  } 
  // genome 0, position 5: background "A"
  {
    variant v{};
    v.position = 5;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome0.add_variant(5, v);
  }

  // genome 1, position 5: background "A"
  {
    variant v{};
    v.position = 5;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome1.add_variant(5, v);
  }
  genome0.process_variant_assignments();
  genome0.get_consensus(0);
  genome1.process_variant_assignments();
  genome1.get_consensus(1);

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

  // SECOND SCENARIO: 5 positions again. Position 2 has two competing alleles,
  // no deletion this time: A@0.90 -> genome 0, T@0.10 -> genome 1, each its own
  // standalone record. Position 3 has a background invariant "A"
  // (half_normal_upper, both genomes) plus a minor insertion +CC@0.10 assigned
  // only to genome 1. All other positions are background "A"s (half_normal_upper,
  // both genomes).
  {
    consensus_sequence genome0(max_position);
    consensus_sequence genome1(max_position);

  // genome 0, position 1: background "A"
  {
    variant v{};
    v.position = 1;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome0.add_variant(1, v);
  }
  // genome 1, position 1: background "A"
  {
    variant v{};
    v.position = 1;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome1.add_variant(1, v);
  }

  // genome 0, position 2 dominant A (freq 0.90)
  {
    variant v{};
    v.position = 2;
    v.nuc = "A";
    v.gapped_freq = 0.90;
    v.assigned_deletion = false;
    v.consensus_numbers = {0};
    genome0.add_variant(2, v);
  }

  // genome 1, position 2 minor variant T (freq 0.10)
  {
    variant v{};
    v.position = 2;
    v.nuc = "T";
    v.gapped_freq = 0.10;
    v.assigned_deletion = false;
    v.consensus_numbers = {1};
    genome1.add_variant(2, v);
  }
  // genome 0, position 3 background "A"
  {
    variant v{};
    v.position = 3;
    v.nuc = "A";
    v.gapped_freq = 1;
    v.assigned_deletion = false;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true; 
    v.consensus_numbers = {0, 1};
    genome0.add_variant(3, v);
  } 
  // genome 1, position 3 background "A"
  {
    variant v{};
    v.position = 3;
    v.nuc = "A";
    v.gapped_freq = 1;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true; 
    v.assigned_deletion = true;
    v.consensus_numbers = {1};
    genome1.add_variant(3, v);
  }
  //genom 1, position 3 minor insertion +CC (freq 0.10)
  {
    variant v{};
    v.position = 3;
    v.nuc = "+CC";
    v.gapped_freq = 0.10;
    v.assigned_deletion = true;
    v.consensus_numbers = {1};
    genome1.add_variant(3, v);
  }
  // genome 0, position 4: background "A"
  {
    variant v{};
    v.position = 4;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome0.add_variant(4, v);
  }
  // genome 1, position 4: background "A"
  {
    variant v{};
    v.position = 4;
    v.nuc = "A";
    v.gapped_freq = 1;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true; 
    v.assigned_deletion = false;
    v.consensus_numbers = {0, 1};
    genome1.add_variant(4, v);
  } 
  // genome 0, position 5: background "A"
  {
    variant v{};
    v.position = 5;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome0.add_variant(5, v);
  }

  // genome 1, position 5: background "A"
  {
    variant v{};
    v.position = 5;
    v.nuc = "A";
    v.gapped_freq = 1.0;
    v.half_normal_upper = true;
    v.position_half_normal_upper = true;
    v.consensus_numbers = {0, 1};
    genome1.add_variant(5, v);
  }

    genome0.process_variant_assignments();
    genome0.get_consensus(2);
    genome1.process_variant_assignments();
    genome1.get_consensus(3);

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

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
