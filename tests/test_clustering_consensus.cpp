//misc. unit tests for clustering code
#include <iostream>
#include <vector>
#include <cmath>
#include "../src/gmm.h"
#include "../src/gmm_1d.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

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

  // FOURTH SCENARIO (TESTs 17, 18, 19, 20 and 21) - gmm_1d::predict prunes a gaussian component holding
  // fewer than -P points (or, when -P is 0, less than -F of the data), moving its points to their
  // next best component. The half normals are never pruned. One fit, reused for every test
  {
    //0.5 is the tiny cluster to prune; 0.04 sits in the left half normal's support
    std::vector<double> x;
    for(int i=0; i < 20; i++) x.push_back(0.285 + i * 0.03 / 19);  //indices 0-19, evenly over 0.285-0.315
    for(int i=0; i < 2; i++) x.push_back(0.5);   //indices 20-21
    for(int i=0; i < 20; i++) x.push_back(0.685 + i * 0.03 / 19);  //indices 22-41, evenly over 0.685-0.715
    for(int i=0; i < 1; i++) x.push_back(0.04);  //indices 42-43
    const uint32_t i03 = 0, i05 = 20, i07 = 22, i004 = 42;

    uint32_t n = 5; //4 gaussians plus the two half normals
    gmm_1d model(n, 42);
    model.set_use_half_normal_for_noise(true, 0.97);
    model.set_covariance_prior(1e-4);
    model.set_mean_precision_prior(1e-2);
    model.set_half_normal_covariance_prior(1e-4);
    model.fit(x);
    const int left_half_normal = n - 2;

    //TEST 17 - the fit places a gaussian at each of 0.3, 0.5 and 0.7
    {
      std::vector<double> means = model.get_means();
      auto has_mean = [&](double target){
        for(double m : means){
          if(std::abs(m - target) < 1e-3) return true;
        }
        return false;
      };
      bool pass = true;
      pass &= check_flag(4, "mean 0.3", "fitted", has_mean(0.3), true);
      pass &= check_flag(4, "mean 0.5", "fitted", has_mean(0.5), true);
      pass &= check_flag(4, "mean 0.7", "fitted", has_mean(0.7), true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 18 - with no pruning, the 0.5 points have their own label and the 0.04 points are in the
    //left half normal
    {
      model.set_min_cluster_points(0);
      model.set_min_cluster_fraction(0.0);
      std::vector<int> labels = model.predict(x);
      bool pass = true;
      pass &= check_flag(4, "0.5", "own label", labels[i05] != labels[i03] && labels[i05] != labels[i07], true);
      pass &= check_flag(4, "0.04", "in left half normal", labels[i004] == left_half_normal, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 19 - -P 3 prunes the 0.5 cluster into a big group, the big groups stay apart, and the left
    //half normal survives with fewer than 3 points
    {
      model.set_min_cluster_points(3);
      model.set_min_cluster_fraction(0.0);
      std::vector<int> labels = model.predict(x);
      bool pass = true;
      pass &= check_flag(4, "0.5", "joined a big group", labels[i05] == labels[i03] || labels[i05] == labels[i07], true);
      pass &= check_flag(4, "0.3 and 0.7", "distinct", labels[i03] != labels[i07], true);
      pass &= check_flag(4, "0.04", "in left half normal", labels[i004] == left_half_normal, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 20 - -F 0.10 prunes the same way when -P is 0 (2 of 44 points is under 10%, 20 is not)
    {
      model.set_min_cluster_points(0);
      model.set_min_cluster_fraction(0.10);
      std::vector<int> labels = model.predict(x);
      bool pass = true;
      pass &= check_flag(4, "0.5", "joined a big group", labels[i05] == labels[i03] || labels[i05] == labels[i07], true);
      pass &= check_flag(4, "0.3 and 0.7", "distinct", labels[i03] != labels[i07], true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 21 - -P overrides -F: -F 0.5 alone would prune both big groups (20 of 44 is under 50%),
    //but with -P 3 only the 0.5 cluster is pruned
    {
      model.set_min_cluster_points(3);
      model.set_min_cluster_fraction(0.5);
      std::vector<int> labels = model.predict(x);
      bool pass = true;
      pass &= check_flag(4, "0.5", "joined a big group", labels[i05] == labels[i03] || labels[i05] == labels[i07], true);
      pass &= check_flag(4, "0.3 and 0.7", "distinct", labels[i03] != labels[i07], true);
      num_tests++;
      if (pass) success++;
    }
  }

  // FIFTH SCENARIO (TESTs 22, 23, 24, 25 and 26) - flag_freq_outlier_variants marks a variant
  // freq_outlier when it sits more than -W standard deviations from the mean of the cluster it was
  // assigned to. Out of range variants and variants with no valid cluster are skipped
  {
    std::vector<double> eff_means = {0.3, 0.7};
    std::vector<double> eff_vars = {0.0004, 0.0004}; //standard deviation 0.02, so -W 3 is 0.06

    variant inside;
    inside.gapped_freq = 0.33; //0.03 from 0.3
    inside.cluster_assigned = 0;

    variant outside;
    outside.gapped_freq = 0.40; //0.10 from 0.3
    outside.cluster_assigned = 0;

    variant out_of_range;
    out_of_range.gapped_freq = 0.99; //far from 0.3, but outside the model's range
    out_of_range.cluster_assigned = 0;
    out_of_range.outside_freq_range = true;

    variant no_cluster;
    no_cluster.gapped_freq = 0.50; //far from both clusters, but assigned to none
    no_cluster.cluster_assigned = -1;

    std::vector<variant> outlier_variants = {inside, outside, out_of_range, no_cluster};
    flag_freq_outlier_variants(outlier_variants, eff_means, eff_vars, 3.0);

    //TEST 22 - a variant inside the cutoff is not an outlier
    {
      bool pass = check_flag(5, "inside", "freq_outlier", outlier_variants[0].freq_outlier, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 23 - a variant past the cutoff is an outlier
    {
      bool pass = check_flag(5, "outside", "freq_outlier", outlier_variants[1].freq_outlier, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 24 - an out of range variant is skipped, however far it is
    {
      bool pass = check_flag(5, "out_of_range", "freq_outlier", outlier_variants[2].freq_outlier, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 25 - a variant with no valid cluster is skipped
    {
      bool pass = check_flag(5, "no_cluster", "freq_outlier", outlier_variants[3].freq_outlier, false);
      num_tests++;
      if (pass) success++;
    }

    //TEST 26 - a variant exactly at the cutoff is not an outlier, the rule is strictly greater than.
    //Run on its own with values binary represents exactly, so the distance lands on the cutoff
    {
      std::vector<double> exact_means = {0.25};
      std::vector<double> exact_vars = {0.0625}; //standard deviation 0.25

      variant at_cutoff;
      at_cutoff.gapped_freq = 0.5; //exactly 0.25 from the mean
      at_cutoff.cluster_assigned = 0;

      std::vector<variant> exact_variants = {at_cutoff};
      flag_freq_outlier_variants(exact_variants, exact_means, exact_vars, 1.0);

      bool pass = check_flag(5, "at_cutoff", "freq_outlier", exact_variants[0].freq_outlier, false);
      num_tests++;
      if (pass) success++;
    }
  }

  // SIXTH SCENARIO (TESTs 27 and 28) - variant_assigner::assign puts an upper half normal variant in
  // every genome and a lower half normal variant in none, whatever cluster or posteriors it carries
  {
    std::vector<double> solution = {0.2, 0.8};
    std::vector<double> means = {0.2, 0.8}; //cluster 0 is genome 0, cluster 1 is genome 1

    variant v_upper;
    v_upper.half_normal_upper = true;
    v_upper.cluster_assigned = 0;
    v_upper.probabilities = {0.9, 0.1}; //would reach genome 0 only without the flag

    variant v_lower;
    v_lower.half_normal_lower = true;
    v_lower.cluster_assigned = 1;
    v_lower.probabilities = {0.1, 0.9}; //would reach genome 1 only without the flag

    variant v_normal;
    v_normal.cluster_assigned = 0;
    v_normal.probabilities = {0.9, 0.1}; //control, routed by its posteriors

    std::vector<variant> assign_variants = {v_upper, v_lower, v_normal};
    variant_assigner(solution, means, 2.0).assign(assign_variants);

    //TEST 27 - the upper half normal variant reaches every genome, while the control with the same
    //posteriors reaches genome 0 only
    {
      bool pass = true;
      pass &= check_flag(6, "v_upper", "in every genome",
                         assign_variants[0].consensus_numbers == std::vector<uint32_t>{0, 1}, true);
      pass &= check_flag(6, "v_upper", "no ambiguous genomes", assign_variants[0].ambiguous_numbers.empty(), true);
      pass &= check_flag(6, "v_normal", "in genome 0 only",
                         assign_variants[2].consensus_numbers == std::vector<uint32_t>{0}, true);
      num_tests++;
      if (pass) success++;
    }

    //TEST 28 - the lower half normal variant reaches no genome
    {
      bool pass = true;
      pass &= check_flag(6, "v_lower", "no confirmed genomes", assign_variants[1].consensus_numbers.empty(), true);
      pass &= check_flag(6, "v_lower", "no ambiguous genomes", assign_variants[1].ambiguous_numbers.empty(), true);
      num_tests++;
      if (pass) success++;
    }
  }

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
