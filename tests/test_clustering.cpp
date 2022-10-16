#include<iostream>
#include "ap.h"
#include "stdafx.h"
#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/clustering.h"
#include "../src/kmeans.h"

void reset_positions(std::vector<position> &all_positions){
  /*
   * @param all_positions : vector containing position objects for whole genome
   * Reset function to basic alleles with no depths for the position vector.
   * Meant to allow multiple tests within a single function.
   */
  all_positions.clear();

  //set up basic nucleotides
  std::vector<std::string> basic_nts = {"A", "C", "G", "T"};
  std::vector<allele> basic_alleles;
  for(std::string nt : basic_nts){
    allele new_allele;
    new_allele.nuc = nt;
    new_allele.depth = 0;
    basic_alleles.push_back(new_allele);
  }
  //populate all positions
  for(uint32_t i=0; i < 29903; i++){
    position new_position;
    new_position.pos = i;
    new_position.depth = 0;
    new_position.ad = basic_alleles;
    all_positions.push_back(new_position);
  }
}
int _integreation_test_consensus(){
  /*
   * Tests the consensus call using two approaches. First, we call consensus using autothresholding.
   * Second, we call consensus using the traditional mpileup and set the same threshold.
   * consensus sequences should be identical.
   */

  return(0);
}

int _integreation_create_haplotypes(){
  return(0);
}

int _unit_test_parse_md_tag(){
  /*
   * @return success : 0 for passing else failing 
   *
   * Tests the process of parsing the MD tag.
   */
  int success = 0;

  /*uint32_t abs_start_pos = 23620;
  uint32_t abs_end_pos = 23727;
  std::vector<position> all_positions;
  std::vector<double> return_frequencies;
  bool reverse = false;

  reset_positions(all_positions);*/

  /* Generate some fake haplotypes / positions to test.
   * Using the amplicon between 23601-23619 23728-23751.
   *
   * Case A. Substituion in the primer region.
   *    Return should be an empty haplotype vector.
   * Case B. Deletion 
   */
  /*uint8_t aux = (uint8_t) 'Z2A7';
  uint8_t seq = (uint8_t) 'ZTTT';
  uint32_t correction_factor = 10;
  uint32_t length = 10;
  std::vector<int> haplotypes;
  std::vector<uint32_t> positions;
  std::vector<uint32_t> ignore_positions;*/
  

  //parse_md_tag(aux, haplotypes, positions, abs_start_pos, all_positions, seq, length, correction_factor, abs_end_pos, ignore_positions, reverse);
  

  return(success);
}

int _unit_test_frequencies(){
  /*
   * Unit test for the function create_frequency_matrix.
   */
  int success = 0;

  int32_t primer_offset = 0;
  std::string pair_info = "../data/contamination_tests/primer_pairs.tsv";
  std::string bed = "../data/contamination_tests/sars_primers_strand.bed";

  std::vector<double> return_frequencies;
  std::vector<primer> primers;
  IntervalTree amplicons;
  std::vector<position> all_positions;
  
  //populate amplicons for sars-cov-2
  primers = populate_from_file(bed, primer_offset);
  amplicons = populate_amplicons(pair_info, primers);
  reset_positions(all_positions);

  /* Generate some fake haplotypes / positions to test.
   * Using the amplicon between 23601-23619 23728-23751.
   *
   * In this fictional genome 23700 reference is T, 23701 reference is also T.
   *
   * Case A. Testing the retention of mutations below/above 0.03 total depth freq.
   *    Drop position 23700 w/ mutation frequency < 0.03, return empty frequency vector
   *    Keep position 23701 w/ mutation frequency @ 0.1, retrun frequency vector w/ 0.10
   * Case B. Relative proportions of mutations in the same position.
   *    Two mutations at 23700 in proportions 0.30, 0.20, match ref 0.50
   * Case C. Handling soft clipped reads.
   *    One mutation at 0.20, 0.60 match ref, and 0.20 soft clipped. SC get eliminated
   *    from the total making final proportions 0.25 mutation 0.75 match ref
   * Case D. Testing the use of a deletion in creation of frequencies.
   *    One deletion at 23700 0.20 with the other 0.80 matching the reference.
   */

  //declare some baseline variables
  std::vector<uint32_t> positions;
  std::vector<int> haplotypes;
  std::vector<double> frequencies;
  std::vector<float> qualities;
  bool reverse = false;
  uint32_t abs_start_pos = 23620;
  uint32_t abs_end_pos = 23727;
  std::vector<std::string> nucleotides;
  std::vector<uint32_t> range = {23620, 23727};
  std::vector<double> frequency_gt = {0.1};
  int i = 0;

  //Case A - Don't return mutations <= 0.03, do return mutations > 0.03 frequency 
  positions = {23700};
  haplotypes = {-2};
  nucleotides = {"T"};
  qualities = {37};
  //add haplotypes that match the reference 100 times
  while(i < 97){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  nucleotides.clear();
  haplotypes.clear();
  haplotypes = {0};
  nucleotides = {"A"};
  i = 0;
  //add a mutation 3 times for 3% frequency
  while(i < 3){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  //add a mutation at over 3% frequency
  positions.clear();
  haplotypes.clear();
  nucleotides.clear();
  positions = {23701};
  nucleotides = {"T"};
  haplotypes = {-2};
  i = 0;
  while(i < 90){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  nucleotides = {"A"};
  haplotypes = {0};
  i = 0;
  //indicates a mutation at 10% frequency
  while(i < 10){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }

  //try to create the frequency matrix
  frequencies = create_frequency_matrix(amplicons, all_positions);
  if(frequencies == frequency_gt){
    success += 1;
  }

  //Case B - Testing mulitple mutations at one position.
  frequency_gt = {0.2, 0.3};
  amplicons.clear();
  reset_positions(all_positions);
  positions = {23700};
  haplotypes = {-2};
  nucleotides = {"T"};
  i = 0;
  //adding 50 match ref, 20 substitution to C, 30 substitution to A
  while(i < 50){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  i = 0;
  haplotypes = {1};
  nucleotides = {"C"};
  while(i < 20){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  i = 0;
  haplotypes = {0};
  nucleotides = {"A"};
  while(i < 30){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  frequencies = create_frequency_matrix(amplicons, all_positions);
  if(frequencies == frequency_gt){
    success += 1;
  }

  //Case C - testing the removal of SC from final read count
  frequency_gt = {0.25};
  amplicons.clear();
  reset_positions(all_positions);
  haplotypes = {-2};
  nucleotides = {"T"};
  i = 0;
  //adding 60 match ref, 20 substitution to C, 20 soft clipped
  while(i < 60){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  i = 0;
  haplotypes = {1};
  nucleotides = {"C"};
  while(i < 20){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  i = 0;
  haplotypes = {-1};
  //don't need to update allele depths for SC
  while(i < 20){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    i++;
  }
  frequencies = create_frequency_matrix(amplicons, all_positions);
  if(frequencies == frequency_gt){
    success += 1;
  }

  //Case D - Test for the use of a deletion
  frequency_gt = {0.20};
  amplicons.clear();
  reset_positions(all_positions);
  haplotypes = {-2};
  nucleotides = {"T"};
  i = 0;
  //adding 80 match ref, 20 deletion
  while(i < 80){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  i = 0;
  haplotypes = {4};
  nucleotides = {"*"}; //this indicate deletion
  while(i < 20){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions, qualities);
    i++;
  }
  frequencies = create_frequency_matrix(amplicons, all_positions);
  if(frequencies == frequency_gt){
    success += 1;
  }
  return(success);
}

int _unit_test_sil_score(){
  /*
   * Makes sure the scoring metric is working to par.
   */
  return(0);  
}

//test the clustering function
int _unit_test_kmeans(){
  /*
   * @return success : 0 for passing else failing
   *
   * Does two checks on the kmeans function looking for (1) accurate cluster centers 
   * and (2) an accurate sil score.
   */

  int success = 0;
  //test case one
  alglib::real_2d_array xy = "[[1],[1],[1],[4],[4],[4]]";
  std::vector<float> ground_truth_centers = {1,4};
  int number_clusters = 2;

  cluster cluster_results;
  k_means(number_clusters, xy, cluster_results);
  //print_cluster_info(cluster_results);
  if(cluster_results.sil_score == 1){
    success += 1;
  }
  //if(cluster_results.centers == ground_truth_centers){
  //  success += 1;
  //}
  //std::cout << "success " << success << std::endl;
  return(success);
}

int main(){
  int success = 4; //count total number of tests
  /*success -= _unit_test_frequencies(); //contains 4 tests
  success -= _unit_test_parse_md_tag(); //contains 2 tests
  success -= _integreation_create_haplotypes();
  
  std::cout << "success: " << success << std::endl;
  */
  //success -= _unit_test_kmeans();
  //success -= _unit_test_sil_score();
  success -= determine_threshold("../data/contamination_tests/test.calmd.bam", "../data/contamination_tests/sequence.fasta", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0, 0.8, 20, 'N', 10, true, "amplicon");
  //determine_threshold("../data/contamination_tests/test.calmd.bam", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0);
  //determine_threshold("../data/contamination_tests/simulated_alpha_beta_90_10.bam", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0);
  return 0;
}
