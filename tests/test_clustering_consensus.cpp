//misc. unit tests for clustering code
#include <iostream>
#include <vector>
#include "../src/gmm.h"
#include "../src/call_consensus_clustering.h"

//prints mismatches in the same format as test_unit_new_consensus
bool check_flag(uint32_t scenario, std::string variant_name, std::string flag_name, bool actual, bool expected){
  if(actual != expected){
    std::cerr << "scenario " << scenario << " " << variant_name << " " << flag_name << " expected " << expected
              << " got " << actual << std::endl;
    return false;
  }
  return true;
}

int main() {
  int num_tests = 0;
  int success = 0;

  // FIRST SCENARIO (TESTs 1, 2, 3 and 4) - set_deletion_flags marks a deletion overlapping a more
  // abundant one as overlapped_deletion, so it stays out of the fit. Deletions are accepted most
  // abundant first; depth flagged deletions take no part
  {
    variant dominant_del;
    dominant_del.position = 10;
    dominant_del.nuc = "-AAA"; //spans 10-12
    dominant_del.gapped_freq = 0.80;

    variant same_start_del;
    same_start_del.position = 10;
    same_start_del.nuc = "-AA"; //spans 10-11, same start as the dominant deletion
    same_start_del.gapped_freq = 0.15;

    variant inside_del;
    inside_del.position = 11;
    inside_del.nuc = "-A"; //starts inside the dominant deletion
    inside_del.gapped_freq = 0.05;

    variant adjacent_del;
    adjacent_del.position = 13;
    adjacent_del.nuc = "-CC"; //spans 13-14, starts right after the dominant deletion ends
    adjacent_del.gapped_freq = 0.40;

    variant low_depth_del;
    low_depth_del.position = 20;
    low_depth_del.nuc = "-GGGG"; //spans 20-23, the most abundant here but depth flagged
    low_depth_del.gapped_freq = 0.90;
    low_depth_del.depth_flag = true;

    variant under_low_depth_del;
    under_low_depth_del.position = 21;
    under_low_depth_del.nuc = "-G"; //inside the depth flagged deletion
    under_low_depth_del.gapped_freq = 0.10;

    std::vector<variant> overlap_variants = {dominant_del, same_start_del, inside_del, adjacent_del,
                                             low_depth_del, under_low_depth_del};
    set_deletion_flags(overlap_variants);

    //TEST 1 - a minor deletion at the same start as the dominant one is overlapped, the dominant is not
    {
      bool pass = true;
      pass &= check_flag(1, "dominant_del", "overlapped_deletion", overlap_variants[0].overlapped_deletion, false);
      pass &= check_flag(1, "same_start_del", "overlapped_deletion", overlap_variants[1].overlapped_deletion, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 2 - a minor deletion starting inside the dominant deletion is overlapped
    {
      bool pass = check_flag(1, "inside_del", "overlapped_deletion", overlap_variants[2].overlapped_deletion, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 3 - a deletion starting right after the dominant one ends does not overlap it
    {
      bool pass = check_flag(1, "adjacent_del", "overlapped_deletion", overlap_variants[3].overlapped_deletion, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 4 - a depth flagged deletion is ignored: it is not flagged, and it does not block a
    //smaller deletion inside it
    {
      bool pass = true;
      pass &= check_flag(1, "low_depth_del", "overlapped_deletion", overlap_variants[4].overlapped_deletion, false);
      pass &= check_flag(1, "under_low_depth_del", "overlapped_deletion", overlap_variants[5].overlapped_deletion, false);
      num_tests++;
      if (pass) success++;
    }
  }

  // SECOND SCENARIO (TESTs 5, 6, 7, 8 and 9) - every base under a deletion writes a reference row
  // at ~1 - the deletion's frequency, so set_deletion_flags keeps only the median reference row of
  // each accepted, in range deletion for the fit and marks the rest deletion_span_duplicate
  {
    //TEST 5 and 9 - an in range deletion over three reference rows, plus a substitution in its span
    variant del10;
    del10.position = 10;
    del10.nuc = "-AAA"; //spans 10-12
    del10.gapped_freq = 0.30;

    variant ref10;
    ref10.position = 10;
    ref10.nuc = "A";
    ref10.ref_allele = true;
    ref10.gapped_freq = 0.70; //median, kept

    variant ref11;
    ref11.position = 11;
    ref11.nuc = "A";
    ref11.ref_allele = true;
    ref11.gapped_freq = 0.69;

    variant ref12;
    ref12.position = 12;
    ref12.nuc = "A";
    ref12.ref_allele = true;
    ref12.gapped_freq = 0.71;

    variant snp11;
    snp11.position = 11;
    snp11.nuc = "T"; //not a reference row
    snp11.gapped_freq = 0.05;

    //TEST 6 - a single base deletion has one reference row under it
    variant del20;
    del20.position = 20;
    del20.nuc = "-C";
    del20.gapped_freq = 0.40;

    variant ref20;
    ref20.position = 20;
    ref20.nuc = "C";
    ref20.ref_allele = true;
    ref20.gapped_freq = 0.60;

    //TEST 7 - depth and quality flagged reference rows take no part in the median
    variant del30;
    del30.position = 30;
    del30.nuc = "-GGGG"; //spans 30-33
    del30.gapped_freq = 0.30;

    variant ref30;
    ref30.position = 30;
    ref30.nuc = "G";
    ref30.ref_allele = true;
    ref30.gapped_freq = 0.70;

    variant ref31;
    ref31.position = 31;
    ref31.nuc = "G";
    ref31.ref_allele = true;
    ref31.gapped_freq = 0.72; //upper of the two eligible rows, kept

    variant ref32;
    ref32.position = 32;
    ref32.nuc = "G";
    ref32.ref_allele = true;
    ref32.gapped_freq = 0.71;
    ref32.depth_flag = true;

    variant ref33;
    ref33.position = 33;
    ref33.nuc = "G";
    ref33.ref_allele = true;
    ref33.gapped_freq = 0.69;
    ref33.qual_flag = true;

    //TEST 8 - a noise level deletion does not depress the reference rows under it
    variant del40;
    del40.position = 40;
    del40.nuc = "-TT"; //spans 40-41
    del40.gapped_freq = 0.02;
    del40.outside_freq_range = true;

    variant ref40;
    ref40.position = 40;
    ref40.nuc = "T";
    ref40.ref_allele = true;
    ref40.gapped_freq = 0.95;

    variant ref41;
    ref41.position = 41;
    ref41.nuc = "T";
    ref41.ref_allele = true;
    ref41.gapped_freq = 0.96;

    std::vector<variant> span_variants = {del10, ref10, ref11, ref12, snp11,
                                          del20, ref20,
                                          del30, ref30, ref31, ref32, ref33,
                                          del40, ref40, ref41};
    set_deletion_flags(span_variants);

    //TEST 5 - the median reference row under the deletion is kept, the other two are duplicates
    {
      bool pass = true;
      pass &= check_flag(2, "ref10", "deletion_span_duplicate", span_variants[1].deletion_span_duplicate, false);
      pass &= check_flag(2, "ref11", "deletion_span_duplicate", span_variants[2].deletion_span_duplicate, true);
      pass &= check_flag(2, "ref12", "deletion_span_duplicate", span_variants[3].deletion_span_duplicate, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 6 - a lone reference row under a deletion is kept
    {
      bool pass = check_flag(2, "ref20", "deletion_span_duplicate", span_variants[6].deletion_span_duplicate, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 7 - the median of the eligible rows is kept, and the flagged rows are never duplicates
    {
      bool pass = true;
      pass &= check_flag(2, "ref30", "deletion_span_duplicate", span_variants[8].deletion_span_duplicate, true);
      pass &= check_flag(2, "ref31", "deletion_span_duplicate", span_variants[9].deletion_span_duplicate, false);
      pass &= check_flag(2, "ref32", "deletion_span_duplicate", span_variants[10].deletion_span_duplicate, false);
      pass &= check_flag(2, "ref33", "deletion_span_duplicate", span_variants[11].deletion_span_duplicate, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 8 - no reference row under a noise level deletion is a duplicate
    {
      bool pass = true;
      pass &= check_flag(2, "ref40", "deletion_span_duplicate", span_variants[13].deletion_span_duplicate, false);
      pass &= check_flag(2, "ref41", "deletion_span_duplicate", span_variants[14].deletion_span_duplicate, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 9 - a non reference variant inside a deletion's span is never a duplicate
    {
      bool pass = check_flag(2, "snp11", "deletion_span_duplicate", span_variants[4].deletion_span_duplicate, false);
      num_tests++;
      if (pass) success++;
    }
  }

  // THIRD SCENARIO (TESTs 10, 11, 12, 13, 14, 15 and 16) - flag_noise_explained_variants moves an
  // allele to the upper half normal when its frequency plus the noise at its position reaches the
  // invariant threshold. Noise is any noise level allele recorded at the position, plus noise level
  // deletions spanning it; insertions and depth flagged rows are not noise
  {
    //TEST 10 - substitution noise
    variant a10;
    a10.position = 10;
    a10.nuc = "A";
    a10.gapped_freq = 0.96;

    variant c10;
    c10.position = 10;
    c10.nuc = "C";
    c10.gapped_freq = 0.02;
    c10.outside_freq_range = true;

    variant g10;
    g10.position = 10;
    g10.nuc = "G";
    g10.gapped_freq = 0.02;
    g10.outside_freq_range = true;

    //TEST 11 - near miss
    variant a20;
    a20.position = 20;
    a20.nuc = "A";
    a20.gapped_freq = 0.90;

    variant c20;
    c20.position = 20;
    c20.nuc = "C";
    c20.gapped_freq = 0.02;
    c20.outside_freq_range = true;

    //TEST 12 - a noise deletion counts at its first base
    variant a30;
    a30.position = 30;
    a30.nuc = "A";
    a30.gapped_freq = 0.96;

    variant del30;
    del30.position = 30;
    del30.nuc = "-TT"; //spans 30-31
    del30.gapped_freq = 0.02;
    del30.outside_freq_range = true;

    //TEST 13 - a noise deletion counts at its later bases, and not past its end
    variant del40;
    del40.position = 40;
    del40.nuc = "-GGG"; //spans 40-42
    del40.gapped_freq = 0.02;
    del40.outside_freq_range = true;

    variant ref41;
    ref41.position = 41;
    ref41.nuc = "G";
    ref41.ref_allele = true;
    ref41.gapped_freq = 0.96;

    variant ref43;
    ref43.position = 43;
    ref43.nuc = "C";
    ref43.ref_allele = true;
    ref43.gapped_freq = 0.96;

    //TEST 14 - an insertion is not noise
    variant a50;
    a50.position = 50;
    a50.nuc = "A";
    a50.gapped_freq = 0.95;

    variant ins50;
    ins50.position = 50;
    ins50.nuc = "+TT";
    ins50.gapped_freq = 0.02;
    ins50.outside_freq_range = true;

    //TEST 15 - a depth flagged row is not noise
    variant a60;
    a60.position = 60;
    a60.nuc = "A";
    a60.gapped_freq = 0.95;

    variant c60;
    c60.position = 60;
    c60.nuc = "C";
    c60.gapped_freq = 0.02;
    c60.outside_freq_range = true;
    c60.depth_flag = true;

    std::vector<variant> noise_variants = {a10, c10, g10,
                                           a20, c20,
                                           a30, del30,
                                           del40, ref41, ref43,
                                           a50, ins50,
                                           a60, c60};
    flag_noise_explained_variants(noise_variants, 0.97);

    //TEST 10 - 0.96 plus 0.04 of substitution noise reaches the threshold
    {
      bool pass = true;
      pass &= check_flag(3, "a10", "half_normal_upper", noise_variants[0].half_normal_upper, true);
      pass &= check_flag(3, "a10", "outside_freq_range", noise_variants[0].outside_freq_range, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 11 - 0.90 plus 0.02 of noise falls short
    {
      bool pass = true;
      pass &= check_flag(3, "a20", "half_normal_upper", noise_variants[3].half_normal_upper, false);
      pass &= check_flag(3, "a20", "outside_freq_range", noise_variants[3].outside_freq_range, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 12 - the noise deletion counts at its first base
    {
      bool pass = true;
      pass &= check_flag(3, "a30", "half_normal_upper", noise_variants[5].half_normal_upper, true);
      pass &= check_flag(3, "a30", "outside_freq_range", noise_variants[5].outside_freq_range, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 13 - the noise deletion counts inside its span, but not past its end
    {
      bool pass = true;
      pass &= check_flag(3, "ref41", "half_normal_upper", noise_variants[8].half_normal_upper, true);
      pass &= check_flag(3, "ref41", "outside_freq_range", noise_variants[8].outside_freq_range, true);
      pass &= check_flag(3, "ref43", "half_normal_upper", noise_variants[9].half_normal_upper, false);
      pass &= check_flag(3, "ref43", "outside_freq_range", noise_variants[9].outside_freq_range, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 14 - the insertion is not counted as noise
    {
      bool pass = true;
      pass &= check_flag(3, "a50", "half_normal_upper", noise_variants[10].half_normal_upper, false);
      pass &= check_flag(3, "a50", "outside_freq_range", noise_variants[10].outside_freq_range, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 15 - the depth flagged row is not counted as noise
    {
      bool pass = true;
      pass &= check_flag(3, "a60", "half_normal_upper", noise_variants[12].half_normal_upper, false);
      pass &= check_flag(3, "a60", "outside_freq_range", noise_variants[12].outside_freq_range, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 16 - frequency plus noise exactly at the threshold is flagged. Run on its own with
    //values binary represents exactly, so the sum lands on the threshold with no rounding
    {
      variant a_exact;
      a_exact.position = 10;
      a_exact.nuc = "A";
      a_exact.gapped_freq = 0.75;

      variant c_exact;
      c_exact.position = 10;
      c_exact.nuc = "C";
      c_exact.gapped_freq = 0.125; //at the noise ceiling, 1 - 0.875
      c_exact.outside_freq_range = true;

      std::vector<variant> exact_variants = {a_exact, c_exact};
      flag_noise_explained_variants(exact_variants, 0.875);

      bool pass = true;
      pass &= check_flag(3, "a_exact", "half_normal_upper", exact_variants[0].half_normal_upper, true);
      pass &= check_flag(3, "a_exact", "outside_freq_range", exact_variants[0].outside_freq_range, true);
      num_tests++;
      if (pass) success++;
    }
  }

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
