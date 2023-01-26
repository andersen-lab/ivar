#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <vector>
#include <numeric>
#include "./alglib/ap.h"
#include "./alglib/dataanalysis.h"
#include "interval_tree.h"
#include "allele_functions.h"
#include "call_consensus_pileup.h"
#include "get_masked_amplicons.h"
#include "primer_bed.h"
#include "./alglib/stdafx.h"
#include "clustering.h"
#include "ref_seq.h"
#include "kmeans.h"
using namespace alglib;

vector<size_t> sort_indexes(const vector<double> &v) {
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

void print_map(std::unordered_map<uint32_t, std::vector<uint32_t>> const &m){
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": ";
        for(uint32_t y : pair.second){
          std::cout << y << " ";
        }
        std::cout << "}\n";
    }
}


void print_map(std::unordered_map<std::string, int> const &m){
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

void print_map(std::unordered_map<int, std::string> const &m)
{
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}
bool check_nucleotide(std::string nt){
  /*
   * Helper function to determine if a nucelotide is a valid string.
   */
  if(nt == "A"){
    return(true);
  }else if (nt == "C"){
    return(true);
  }else if (nt == "G"){
    return(true);
  }else if(nt == "T"){
    return(true);
  }else if(nt == "N"){
    return(true);
  }else{
    return(false);
  }
}

std::vector<std::vector<uint32_t>> transpose(const std::vector<std::vector<uint32_t>> data) {
  // this assumes that all inner vectors have the same size and
  // allocates space for the complete result in advance
  std::vector<std::vector<uint32_t> > result(data[0].size(),
                                        std::vector<uint32_t>(data.size()));
  for (std::vector<uint32_t>::size_type i = 0; i < data[0].size(); i++) 
      for (std::vector<uint32_t>::size_type j = 0; j < data.size(); j++) {
          result[i][j] = data[j][i];
      }
  return result;
}


std::vector<std::vector<int>> transpose(const std::vector<std::vector<int>> data) {
  // this assumes that all inner vectors have the same size and
  // allocates space for the complete result in advance
  std::vector<std::vector<int> > result(data[0].size(),
                                        std::vector<int>(data.size()));
  for (std::vector<int>::size_type i = 0; i < data[0].size(); i++) 
      for (std::vector<int>::size_type j = 0; j < data.size(); j++) {
          result[i][j] = data[j][i];
      }
  return result;
}

void remove_reference_matches(std::vector<uint32_t> &positions, std::vector<std::vector<int>> &haplotypes){
  /*
   * Takes in a zipped pair of positions and haplotypes and removes
   * the locations which are universally SC or match the ref.
   */
  std::vector<std::vector<int>> transposed_vector = transpose(haplotypes);
  std::vector<std::vector<int>> final_haplotypes;
  std::vector<uint32_t> final_positions;
  for(uint32_t y = 0; y < transposed_vector.size(); y++){
    //check if all values are negative
    bool zeros = std::all_of(transposed_vector[y].begin(), transposed_vector[y].end(), [](int i) { return i< 0; });
    if(!zeros){
      final_haplotypes.push_back(transposed_vector[y]);
      final_positions.push_back(positions[y]);
    }
  }

  positions.clear();
  haplotypes.clear();
  if(final_positions.size() > 0 && final_haplotypes.size() > 0){
    std::vector<std::vector<int>> transposed_haplotypes = transpose(final_haplotypes); 
    positions = final_positions;
    haplotypes = transposed_haplotypes;
  }
}

std::vector<uint32_t> flatten(std::vector<std::vector<uint32_t>> const &vec)
{
    std::vector<uint32_t> flattened;
    for (auto const &v: vec) {
        flattened.insert(flattened.end(), v.begin(), v.end());
    }
    return flattened;
}

void zip(
    const std::vector<int> &haplotypes, 
    const std::vector<uint32_t> &positions, 
    std::vector<std::pair<uint32_t,int>> &zipped)
{
    for(size_t i=0; i< positions.size(); ++i)
    {
        zipped.push_back(std::make_pair(positions[i], haplotypes[i]));
    }
}

void unzip(const std::vector<std::pair<uint32_t, int>> &zipped, std::vector<int> &haplotypes, std::vector<uint32_t> &positions){
  for(size_t i=0; i < positions.size(); i++){
    positions[i] = zipped[i].first;
    haplotypes[i] = zipped[i].second;
  }
  positions.resize(zipped.size());
  haplotypes.resize(zipped.size());
}

void reorder_haplotypes(std::vector<int> &haplotypes, std::vector<uint32_t> &positions){
  /*
   * @param haplotypes : the hapltyope encoded nts
   * @param positions : the position relative to the reference
   *
   * Function reorders both vectors so that the positions are in relative order.
   */

    // Zip the vectors together
    std::vector<std::pair<std::uint32_t, int>> zipped;
    zip(haplotypes, positions, zipped);

    // Sort the vector of pairs
    std::sort(std::begin(zipped), std::end(zipped));

    // Write the sorted pairs back to the original vectors
    unzip(zipped, haplotypes, positions);

}

//print the cluster infor
void print_cluster_info(cluster cluster_results){
  /*
   * @param cluster_results : the cluster object containing clustering stats
   */
  std::cout << "N Clusters: " << cluster_results.n_clusters << std::endl;
  std::cout << "Sil Score: " << cluster_results.sil_score << std::endl;
  for(double x: cluster_results.centers){
    std::cout << "center @: " << x << std::endl;
  }
}

std::string decoded_nucs(int tmp, std::unordered_map<int, std::string> dict_decoded){
  /*
   * @param tmp : integer value encoding for a nucleotide
   * @return dencoded_nuc : nucleotide sequence as a string
   *
   * TODO could be a case-switch?
   */
  std::string decoded_nuc = dict_decoded[tmp];
  return(decoded_nuc);
}


int encoded_nucs(std::string &tmp, std::unordered_map<std::string, int> &dict_encode, std::unordered_map<int, std::string> &dict_decode){
  /*
   * @param tmp : nucleotide to encode
   * @return encoded_nuc : nucleotide encoded as an int
   */
  int encoded_nuc=-200;
  std::unordered_map<std::string, int>::iterator it;
  auto map_it = dict_encode.find(tmp);
  if(map_it == dict_encode.end()){
    int max = -3;
    for(it = dict_encode.begin(); it != dict_encode.end(); it++){
      if(it->second > max){
        max = it->second;
      }
    }  
    encoded_nuc = max + 1;
    dict_decode[max+1] = tmp;
    dict_encode[tmp] = max + 1;
  }else{
    encoded_nuc = dict_encode[tmp];
  }
  return(encoded_nuc);
}

void parse_md_tag(uint8_t *aux, std::vector<int> &haplotypes, std::vector<uint32_t> &positions, std::vector<float> &saved_qualities,
    uint32_t abs_start_pos, std::vector<position> &all_positions,
    uint8_t *seq, uint32_t length, uint32_t correction_factor, uint32_t abs_end_pos,
    std::vector<uint32_t> ignore_positions, bool reverse, bam1_t *r, uint8_t *qualities,
    std::unordered_map<std::string,int> &dict_encode, std::unordered_map<int, std::string> &dict_decode,
    std::vector<uint32_t> insertion_start, std::vector<uint32_t> insertion_length, ref_antd reference, std::string region_){
  /*
   * @param aux : the md tag
   * @param haplotypes : vector with encoded nuc haplotypes
   * @param positions : the positions of sub,ins, and dels that make the haplotype
   * @param abs_start_pos : the start position relative to the reference (is actually abs_end_pos) if reverse
   * @param ad : vector containing variant alleles
   * @param seq : pointer to the sequence
   *
   * Parse the MD tag, populating the haplotypes and positions of the amplicon with
   * substitutions.
   * Any discrepency between the cigar string and MD tag means insertions are in the read.
   * Cigar code handles insertions elsewhere.
   * Qualities are indexed from 0 including the SC region
   */
  //also record the pos that match the reference
  std::vector<uint32_t> ref_pos;
  std::vector<int> ref_haplotypes;
  std::vector<float> ref_qual;

  bool deletion = false; //helping track the deletions
  std::string digits; //helping to track number operations 
  std::string nucs;
  std::string nt;
  std::vector<uint32_t> deletions; //record deletion spots to skip when adding ref nucs
  uint32_t del_correction_factor = 0;
  std::vector<uint32_t> rinsertion_length = insertion_length;
  std::vector<uint32_t> rinsertion_start = insertion_start;
  //int qual = 0;
  int i = 0;
  std::vector<std::string> nucleotides; //store the substitutions & deletions

  do {
    char tmp = aux[i]; //this is the reference nuc 
    //add insertion corrections as needed
    uint32_t insertion_p = 0;
    uint32_t insertion_l = 0;
    for(uint32_t x=0; x < insertion_start.size(); x++){
      if(abs_start_pos >= insertion_start[x]){
        insertion_l = insertion_length[x];
        insertion_p = insertion_start[x];
        abs_start_pos += insertion_l;
        break;
      }
    }
    if(insertion_p != 0 && insertion_l != 0){
      insertion_start.erase(std::remove(insertion_start.begin(), insertion_start.end(), insertion_p), insertion_start.end());
      insertion_length.erase(std::remove(insertion_length.begin(), insertion_length.end(), insertion_l), insertion_length.end());
   }

    if(isdigit(tmp)){ //on digit character
      if(deletion){
        //del_correction_factor += nucs.length(); //this speaks to the relative position
        deletion = false;
        abs_start_pos += stoi(digits);
        for(uint32_t len = 0; len < nucs.length(); len++){
          deletions.push_back(abs_start_pos + len);
          positions.push_back(abs_start_pos + len);
          nt = '*';
          haplotypes.push_back(encoded_nucs(nt, dict_encode, dict_decode));
          nucleotides.push_back(nt);
          saved_qualities.push_back(0);
        }
        abs_start_pos += nucs.length();
        nucs.clear();
        digits.clear();
      }
      deletion = false;
      digits += tmp;
    }else if (isalpha(tmp)){ //on alpha character, not the very first one
      if(deletion){ //we already know this will be a deletion
        nucs += tmp;
        i++;
        continue;
      }
      if(digits != ""){
        if(nucs == ""){
          abs_start_pos += stoi(digits);
          digits.clear();
          continue;
        }
        //check if the position already got soft clipped
        abs_start_pos += std::stoi(digits) + nucs.length();
        nt = "";
        nt = seq_nt16_str[bam_seqi(seq, abs_start_pos+correction_factor - length - 1 - del_correction_factor)];
        char ref = reference.get_base(abs_start_pos, region_);
        if(ref == nt[0]){
          nucs.clear();
          digits.clear();
          continue;
        }
        positions.push_back(abs_start_pos);
        saved_qualities.push_back(qualities[abs_start_pos+correction_factor - length - 1 - del_correction_factor]+0);
        nucleotides.push_back(nt);
        haplotypes.push_back(encoded_nucs(nt, dict_encode, dict_decode));
        digits.clear();
        nucs.clear();
        nucs += tmp;

      }else{
        nucs += tmp;
      }
    }else if (tmp == '^'){
      nucs.clear(); //deletions recorded after
      deletion = true; 
    }
    i++;
  } while(aux[i] != '\0');

  int soft_clipped = -2;
  uint32_t seq_pos = 0;
  std::vector<uint32_t>::iterator it;
  std::vector<uint32_t>::iterator it_ignore;
  std::vector<uint32_t>::iterator it_deletion;

  //fill out the reference matching positions too
  //main issue to watch out for is deletions will shift the NT pos
  uint32_t relative_seq_pos = length;
  uint32_t add_correct = 0;
  for(uint32_t z = length+1; z <= abs_end_pos; z++){
    //add insertion corrections as needed
    uint32_t rinsertion_p = 0;
    uint32_t rinsertion_l = 0;
    for(uint32_t x=0; x < rinsertion_start.size(); x++){
      if(z >= rinsertion_start[x]){
        rinsertion_l = rinsertion_length[x];
        rinsertion_p = rinsertion_start[x];
        relative_seq_pos += rinsertion_l;
        break;
      }
    }
    if(rinsertion_p != 0 && rinsertion_l != 0){
      rinsertion_start.erase(std::remove(rinsertion_start.begin(), rinsertion_start.end(), rinsertion_p), rinsertion_start.end());
      rinsertion_length.erase(std::remove(rinsertion_length.begin(), rinsertion_length.end(), rinsertion_l), rinsertion_length.end());
   }
   //if this is a deletion position we pass it
    it_deletion = std::find(deletions.begin(), deletions.end(), z);
    if(it_deletion != deletions.end()){
      continue;
    }
    it = std::find(positions.begin(), positions.end(), z);
    it_ignore = std::find(ignore_positions.begin(), ignore_positions.end(), z);
    //position wasn't found to be a variant and wasn't found to be deletion
    if(it == positions.end() && it_ignore == ignore_positions.end()){
      ref_pos.push_back(z);
      seq_pos = relative_seq_pos - length + correction_factor + add_correct;
      nt = seq_nt16_str[bam_seqi(seq, seq_pos)];
      if(check_nucleotide(nt) == false){
        std::cout << bam_get_qname(r) << " " << abs_start_pos << " " << abs_end_pos << " " << aux <<" " << nt <<  " " << relative_seq_pos - length + correction_factor + add_correct << std::endl;
      }
      ref_qual.push_back(qualities[seq_pos] + 0);
      int tmp2 = encoded_nucs(nt, dict_encode, dict_decode);
      ref_haplotypes.push_back(tmp2);
      relative_seq_pos += 1;
    }else{
      add_correct++;
    }
  }
  

  if(ref_pos.size() > 0){
    update_allele_depth(all_positions, ref_haplotypes, ref_pos, ref_qual, dict_decode);
  }
  if(positions.size() > 0){
    update_allele_depth(all_positions, haplotypes, positions, saved_qualities, dict_decode);
  }
  for(uint32_t rp : ref_pos){
    positions.push_back(rp);
    haplotypes.push_back(soft_clipped);
  }
}


//calculate the cluster centers
void calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters, cluster &cluster_results){
  /*
   * @params X : array containing all data points
   * @params rep : kmeans reporter object  
   * @returns centers : vector containing all center values 
   */
  std::vector<double> centers;
  double summation [n_clusters][2]; //track all points per cluster
  double point=0;
  int label=0;

  for (int i = 0; i < X.rows(); i++) {
    point = X[i][0];
    label = rep.cidx[i];
    summation[label][0] += point;
    summation[label][1] += 1;
  }
  
  for (int k = 0; k < n_clusters; k++) {
    centers.push_back(summation[k][0]/summation[k][1]); 
  }
  cluster_results.centers = centers;

}

//support function for sil score
double cluster_point_distances(alglib::real_2d_array X, alglib::kmeansreport rep, double point, double center, int n_clusters){
  /*
  * @params X : array containing all data points
  * @params rep : kemans reporter object
  * @params point : the data point of interest
  * @params center : the center of the data point of interests cluster
  * @returns ab : tuple with points a and b representing 
  */
  std::vector<double> internal_dist; //keep track of distance between point of interest and all points in cluster
  double external_dist [n_clusters][2]; //dim 1 cluster value, dim 2 distance sum points and point count
  double compare_point=0;
  int compare_label=0;
  double a=0; //sil score a term
  double b=1; //sil score b term
  double tmp=0;
  //iterate through and mind distances between points in cluster and out of clusters
  for (int i = 0; i < X.rows(); i++) {
    compare_point = X[i][0];
    compare_label = rep.cidx[i];
    //belongs to the same cluster
    if(compare_label == center){
      internal_dist.push_back(abs(point-compare_point));
    }else{ //doesn't belong to the same cluster
      external_dist[compare_label][0] += abs(point-compare_point); //sum distance between poi and cluster points
      external_dist[compare_label][1] += 1; //count number points in cluster    
    }
  } 
  a = average(internal_dist);

  //find the average distance from each external cluster
  for(int z=0; z < n_clusters; z++){
    double zz = z;
    if(zz == center){
      continue;
    }
    tmp = external_dist[z][0] / external_dist[z][1];
    //find the minimum of the external cluster dists
    if(tmp < b){
      b = tmp;
    }
  }
  return((b - a) / std::max(a,b));
}

//function calculates the sil score of all points
void calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep, 
    int n_clusters, cluster &cluster_results){
  /*
  * @params X : array containing data points
  * @params rep : kmeans reported object
  * @params n_clusters : the number of clusters
  * @params cluster_results : the object for storing clustering metrics
  *
  * Calculate the silohuette score for each sample, and then average them
  * together to find the cumulative score.
  * (b - a) / max(a, b) sil score formula
  */
  int rows = X.rows();
  double point = 0;
  double center = 0;
  std::vector<double> sil_scores;
  std::vector<std::vector<double>> sorted_points(n_clusters);
  double tmp = 0;
  for (int i = 0; i < rows; i++) {
    point = X[i][0];
    center = rep.cidx[i];
    sorted_points[center].push_back(point);
    tmp = cluster_point_distances(X, rep, point, center, n_clusters);
    sil_scores.push_back(tmp);
  }

  //store the cumulative results
  cluster_results.sil_score = average(sil_scores);
  cluster_results.sil_scores = sil_scores;
  cluster_results.sorted_points = sorted_points;
}
//assigns per cluster an upper and lower bound
void find_cluster_bounds(cluster &cluster_results){
  /*
   * Function takes each cluster and records the upper and lower bounds. Only needed for upper 
   * cluster technically but recording it for every cluster, for future use.
   */
  std::vector<double> bounds; //temp holding for cluster bounds
  double min;
  double max;
  std::vector<double> tmp;
  std::vector<double> clean_cluster;
  double variance = 0;
  double std_deviation = 0;
  double outlier_lower = 0;
  double outlier_upper = 0;
  double avg = 0;
  for(uint32_t i = 0; i < cluster_results.sorted_points.size(); i++){
    variance = 0;
    std_deviation = 0;
    outlier_lower = 0;
    outlier_upper = 0;
    avg = 0;
    clean_cluster.clear();
    
    tmp = cluster_results.sorted_points[i];
    //look for an remove outliers 
    avg = average(tmp);
    for(uint32_t i = 0; i < tmp.size(); i++){
      //std::cout << "before: " << tmp[i] << "\n";
      variance += std::pow(tmp[i] - avg, 2);
    }
    variance = variance / tmp.size();
    std_deviation = std::sqrt(variance);
    outlier_lower = avg - (std_deviation * 2);
    outlier_upper = avg + (std_deviation * 2);
    //std::cout << "average: " << avg << " std dev: " << std_deviation << std::endl;
    //std::cout << "lower: " << outlier_lower << " upper: " << outlier_upper << std::endl;

    for(uint32_t i = 0; i < tmp.size(); i++){
      if(tmp[i] >= outlier_lower && tmp[i] <= outlier_upper){
        clean_cluster.push_back(tmp[i]);
      }else{
        //clean_cluster.push_back(tmp[i]);
        //std::cout << "outlier: " << tmp[i] << std::endl;
      }
    }
    //std::cout << "cleaned cluster size " << clean_cluster.size() << std::endl;
    /*for(uint32_t i = 0; i < clean_cluster.size(); i++){
      std::cout << clean_cluster[i] << std::endl;
    }*/
    //these find the min/max index within the vector
    min = std::min_element(clean_cluster.begin(), clean_cluster.end()) - clean_cluster.begin();
    max = std::max_element(clean_cluster.begin(), clean_cluster.end()) - clean_cluster.begin();
    bounds = {clean_cluster[min], clean_cluster[max]};
    cluster_results.cluster_bounds.push_back(bounds);
    bounds.clear();
  }
}

//does the actual k means ++ clustering in a loop
void k_means(int n_clusters, alglib::real_2d_array xy, cluster &cluster_results){
  /*
  * @params n_clusters : then number of clusters
  * @params xy : matrix to do kmeans clustering on
  * @params cluster_results : object storing the output of k means clustering
  *
  * Performs K means ++ clustering.
  */
  alglib::clusterizerstate s;
  alglib::kmeansreport rep;
 
  cluster_results.n_clusters = n_clusters;
  int num_points = xy.rows(); //number of points in the data
  clusterizercreate(s);
  //X is the data
  //the number of pointss
  //last pos is distance matrix type
  //num_features is next
  //then distance matrix
  clusterizersetpoints(s, xy, num_points, 1, 2);
  clusterizersetkmeanslimits(s, 5, 0);
  
  //this is the cluster size!!!
  custom_kmeans(s, n_clusters, rep);
  
  if (int(rep.terminationtype) != 1){
    std::cout << "Error in clustering haplotypes" << std::endl;
    exit(1);
  }
  calculate_sil_score(xy, rep, n_clusters, cluster_results);
  calculate_cluster_centers(xy, rep, n_clusters, cluster_results);
  find_cluster_bounds(cluster_results);
}

int count_nucleotides_aux(uint8_t *aux){
  /*
   * @param aux : auxilary tag for alignment
   * @return aux_count : the number of nucleotides account for by the aux tag
   *
   * Given an aux tag, count the number of nucleotides its accounts for.
   */
  int aux_count = 0;
  int i = 0;
  char tmp;
  bool deletion = false;
  std::string digits = "";
  std::string nucs = "";
  do {
    tmp = aux[i]; 
    if(isdigit(tmp)){ //on digit character
      if(deletion){
        deletion = false;
        //std::cout << "add deletion " << digits << std::endl;
        aux_count += stoi(digits);
        nucs.clear();
        digits.clear();
        continue;
      }
      deletion = false;
      digits += tmp;
    }else if(isalpha(tmp)){ //on alpha character, not the very first one
      if(deletion){ //we already know this will be a deletion
        i++;
        continue;
      }
      if(digits != ""){
        if(nucs == ""){
          //std::cout << "add digits a " << stoi(digits) << std::endl;
          aux_count += stoi(digits);
          digits.clear();
          continue;
        }
        if(nucs == "Z"){
          nucs.clear();
        }
        //check if the position already got soft clipped
        //std::cout << "add digits " << std::stoi(digits) + nucs.length() << std::endl;
        aux_count += std::stoi(digits) + nucs.length();
        digits.clear();
        nucs.clear();
        nucs += tmp;
      }else{
        nucs += tmp;
      }
    }else if (tmp == '^'){
      if(nucs != "Z"){
        aux_count += nucs.length();
      }
      //std::cout << "add digits c " << nucs.length() << std::endl; 
      nucs.clear(); //deletions recorded after
      deletion = true; 
    }
    //std::cout << aux[i] << " " << digits << " " << nucs << " " << std::endl;
    
    i++;
  }while(aux[i] != '\0');
  if(nucs != "Z"){
    //std::cout << "add digits e " << stoi(digits) + nucs.length() << std::endl;
    aux_count += stoi(digits) + nucs.length();  
  }else{
    //std::cout << "add digits f " << digits << std::endl;
    aux_count += stoi(digits);
  }
  return(aux_count);
}

void iterate_reads(bam1_t *r, IntervalTree &amplicons, std::vector<position> &all_positions, ref_antd reference, std::string region_, std::unordered_map<std::string, int> &dict_encode, std::unordered_map<int, std::string> &dict_decode){
  /*
   * @param r : alignment object
   * 
   */

  //get the cigar operation 
  uint32_t *cigar = bam_get_cigar(r);
  uint8_t *qualities = bam_get_qual(r);
  uint32_t i = 0;  

  //keep track of operation and length of operation
  uint32_t op = 0;
  uint32_t op_len = 0;

  //we record the range of the read because if thye don't overlap perfectly could
  //accidentally falsely create multiple haplotypes
  std::vector<uint32_t> range;

  //1 = A, 2 = C, 8 = G, 15 = N bits 
  //get a pointer to the sequence
  uint8_t *seq = bam_get_seq(r);

  //get pointer to MD:Z tag which tells us where substitutions are
  uint8_t *aux = bam_aux_get(r, "MD");
  if(aux == NULL){
    std::cout << "No md tag generated." << std::endl;
  }
  //temp variable to count positions
  uint32_t start = 0;
  //this refers to the start position relative to the reference
  bool reverse = bam_is_rev(r);
  uint32_t abs_start_pos = r->core.pos; //leftmost coordinate on ref
  uint32_t abs_end_pos  = bam_endpos(r); //rightmost coordinate on ref

  range.push_back(abs_start_pos); //record the read range
  range.push_back(abs_end_pos);

  //these will later be place in the amplicon object
  std::vector<int> haplotypes;
  std::vector<uint32_t> positions;
  std::vector<float> saved_qualities;
  std::vector<uint32_t> ignore_positions; //all positions that are soft clipped
  std::string nucs = "";
  uint32_t correction_factor = 0;
  bool first_pass = true;
  bool second_pass = true;
  
  //track the length and start of all insertions in the sequence
  std::vector<uint32_t> insertion_start;
  std::vector<uint32_t> insertion_length;

  char nt = 0;
  char ref = 0; //reference base at this pos
  bool primer_mutation = false; //track whether this read has a primer mut
  std::string qname = bam_get_qname(r);
  std::string test = "A01535:8:HJ3YYDSX2:4:1502:9037:35806";
  i = 0;
  int total_nucleotides = 0; //check to make sure the aux matches the cigar
  //iterate through cigar ops for this read
  while(i < r->core.n_cigar){   
    op = bam_cigar_op(cigar[i]); //cigar operation
    op_len = bam_cigar_oplen(cigar[i]); //cigar length
    /*if(qname == test){
      std::cout << op << " " << op_len << " " << aux << std::endl;
    }*/
    //consumes the query but not soft clipped
    if((bam_cigar_type(op) & 1) && (op != 4) && (op != 1)){
      total_nucleotides += op_len;
    }    
    if(op == 4 && first_pass){
      if(!reverse){
        correction_factor = op_len;
      }else{
        correction_factor = op_len;
      }
      first_pass = false;
      ignore_positions.push_back(abs_start_pos+start);
    }
    //figure out if this is the SC region associated with trimming, if so look for mutations
    if((op == 4 && second_pass && !reverse) || (op == 4 && reverse && i+1 == r->core.n_cigar)){
      second_pass = false;     
      uint32_t starting_pos = 0;
      if(reverse){starting_pos += correction_factor;}
      //look at either side regardless of clipped or not clipped
      for(uint32_t pos = start+starting_pos; pos < start+op_len+starting_pos; pos++){
        //if the query begins before the reference we ignore it
        if(pos+abs_start_pos < correction_factor+1){
          continue;
        }
        nt = seq_nt16_str[bam_seqi(seq, pos)];
        ref = reference.get_base(pos+abs_start_pos-correction_factor+1, region_);
        if(nt != ref){
          primer_mutation = true;
        }
      }
    }
    if(op != 4){ first_pass = false;}
    //these are the only operations we care about
    if(op == 1){
      nucs.clear();
      uint32_t start_insertion = start+correction_factor;
      float qual = 0;
      //go get each nt in the insertion region
      for(uint32_t x = 0; x < op_len; x++){
        //the nucelotide at the insertion poi
        nt = seq_nt16_str[bam_seqi(seq, start_insertion+x)];
        if(x == 0){
          nucs += '+';
        }
        nucs += nt;
        qual += qualities[start_insertion+x];
       }
       insertion_start.push_back(abs_start_pos+start);
       insertion_length.push_back(nucs.size()-1);
       haplotypes.push_back(encoded_nucs(nucs, dict_encode, dict_decode));
       positions.push_back(abs_start_pos+start);
       saved_qualities.push_back(qual/nucs.size());
    }

    if (bam_cigar_type(op) & 2){
      start += op_len; //total positions iterated relation to ref
    }                   
    i++;
  }
  //if(qname == test){
  int aux_nt_count = count_nucleotides_aux(aux);
  if(aux_nt_count != total_nucleotides){
    i = 0;
    while(i < r->core.n_cigar){   
      op = bam_cigar_op(cigar[i]); //cigar operation
      op_len = bam_cigar_oplen(cigar[i]); //cigar length
      std::cout << op << " " << op_len << std::endl;                                          
      i++;                                          
    }
    return;
  }
  //}
  parse_md_tag(aux, haplotypes, positions, saved_qualities, abs_start_pos, all_positions, seq, abs_start_pos, correction_factor, abs_end_pos, ignore_positions, reverse, r, qualities, dict_encode, dict_decode, insertion_start, insertion_length, reference, region_);
  if(positions.size() > 0){
    //reoder the positions to be consistent
    reorder_haplotypes(haplotypes, positions);
    
    //places haplotype on amplicon node
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions, primer_mutation, r); 
  }
}

void count_haplotype_occurences(std::vector<std::vector<int>> all_haplotypes, std::vector<std::vector<int>> &save_haplotypes, std::vector<double> &save_read_counts, double &adjusted_read_count, std::vector<uint32_t> final_positions, std::vector<double> &subtract_read_counts){
  /*
   * Combine haplotypes that are the same, count the unique ones, ajust the read counts to remove things that only include soft clipping and matching the reference.
   */
  //find the unique haplotypes in the transposed set
  std::vector<std::vector<int>> unique_haplotypes;
  std::vector<double> count_haplotypes;
  std::vector<uint32_t> flat_pairs;

  for(std::vector<int> exp_haplo : all_haplotypes){
    std::vector<std::vector<int>>::iterator it = std::find(unique_haplotypes.begin(), unique_haplotypes.end(), exp_haplo);
    bool zeros = std::all_of(exp_haplo.begin(), exp_haplo.end(), [](int i) { return i < 0; }); 
    bool mr = std::all_of(exp_haplo.begin(), exp_haplo.end(), [](int i) { return i == -2; }); 
   /*if (std::count(exp_haplo.begin(), exp_haplo.end(), -1)){
      continue;
    }*/
    if(zeros && !mr){
     continue;
    }
    //adjusted_read_count += 1;
    //found
    if(it != unique_haplotypes.end()){
      int index = it - unique_haplotypes.begin();
      count_haplotypes[index] += 1;
      subtract_read_counts[index] += 1;
    }else{
      unique_haplotypes.push_back(exp_haplo);
      count_haplotypes.push_back(1);
      subtract_read_counts.push_back(1);
    }
  }
  

  final_positions[0] = final_positions[0];
  bool match = true;
  int position_1;
  int position_2;
  std::vector<uint32_t>::iterator it_3;
  std::vector<uint32_t>::iterator it_2;
  std::vector<int>::iterator it;
  std::vector<std::vector<uint32_t>> pairs;
  std::unordered_map<uint32_t, std::vector<uint32_t>> matches;
  //here we look to combine haplotypes with soft-clipped or not covered regions (-1) with existing haplotypes
  for(uint32_t i = 0; i < unique_haplotypes.size(); i++){
    if(count_haplotypes[i] < 10){continue;}
    for(uint32_t c = 0; c < unique_haplotypes.size(); c++){
     if(count_haplotypes[c] < 10){continue;}
     if(i == c){continue;} //identity
      //if this isn't doft clipped it can't be a match
      it = std::find(unique_haplotypes[c].begin(), unique_haplotypes[c].end(), -1);
      if(it == unique_haplotypes[c].end()){
        continue;
      }
            
      match = true; 
      //co-iterate the two haplotypes to compare breaking as soon as we hit a differing NT != -1
      for(uint32_t z = 0; z < unique_haplotypes[i].size(); z++){
        position_1 = unique_haplotypes[i][z];
        position_2 = unique_haplotypes[c][z];        
        //we differ at a non soft-clipped or not covered location
        if(position_1 != -1 && position_2 != -1 && position_1 != position_2){
          match = false;
        }
        if(!match){
          break;
        }
      }
      if(match){
        matches[i].push_back(c);
      }
    }
  }
  
  std::vector<uint32_t> tmp;
  for(auto i: sort_indexes(count_haplotypes)) {
    tmp.clear();
    tmp.push_back(i);
    flat_pairs.push_back(i);
    for(uint32_t y : matches[i]){
      //it = std::find(flat_pairs.begin(), flat_pairs.end(), y);
      //if(it ==flat_pairs.end()){
        tmp.push_back(y);
        flat_pairs.push_back(y);
      //}    
    }
    pairs.push_back(tmp);  
  }

  //the first value is the softclipped haplotype, the second is the match 
  for(std::vector<uint32_t> p : pairs){
    save_haplotypes.push_back(unique_haplotypes[p[0]]);
    save_read_counts.push_back(count_haplotypes[p[0]]);
    adjusted_read_count += count_haplotypes[p[0]];
    for(uint32_t y=0; y < p.size(); y++){
      if(y == 0){continue;}
      save_read_counts[save_read_counts.size()-1] += count_haplotypes[p[y]];
      adjusted_read_count += count_haplotypes[p[y]];
    }
   }

}

void check_primer_binding_issues(std::vector<uint32_t> final_positions, std::vector<position> &all_positions, std::vector<double> save_read_counts, double adjusted_read_count, std::vector<std::vector<int>> final_haplotypes, std::vector<uint32_t> &suspect_positions, std::vector<uint32_t> &masked_positions, std::vector<int> &masked_alleles, std::unordered_map<int, std::string> dict_decode, bool primer_issue, bool &flagged, std::vector<double> subtract_read_counts, ref_antd reference, std::string region_, std::vector<double> &frequencies){
  /*
   * If we know this primer is masked, we do an additional check in which we make sure the mutation
   * frequencies for all haplotypes are close to the total depth frequencies.
   */
  suspect_positions.clear();

  std::vector<std::vector<int>> transposed_haplotypes = transpose(final_haplotypes);
 
  std::vector<int> mut; //tracks the mutations
  std::vector<double> mut_count; //tracks how many times we've seen the mutation
  std::vector<double> sub_mut_count;                                 
  std::vector<std::vector<uint32_t>> index; //track index of haplotypes where mutation is found
  std::vector<int>::iterator it;
  
  double arc = 0;
  for(uint32_t i = 0; i < final_positions.size(); i++){
    mut.clear();
    mut_count.clear();
    index.clear();
    arc = 0;
    //look for this position to not match ref in haplotypes
    std::vector<int> all_versions = transposed_haplotypes[i];
    
    //this iterates through the different values found at this position for all haplotypes
    char ref = reference.get_base(final_positions[i], region_);
    std::vector<int> tmp_all_versions;
    for(int i : all_versions){
      if(i != -2){
        tmp_all_versions.push_back(i);
      }else{
        if(ref == 'A'){
          tmp_all_versions.push_back(0);
        }else if(ref == 'C'){
          tmp_all_versions.push_back(1);
        }else if(ref == 'G'){
          tmp_all_versions.push_back(2);
        }else if(ref == 'T'){
          tmp_all_versions.push_back(3);
        }
      }
     }
    for(uint32_t x = 0; x < all_versions.size(); x++){
      if(all_versions[x] == -2 || all_versions[x] >= 0){
        arc += save_read_counts[x];
      }
      //don't apply to deletions
      if((tmp_all_versions[x] >= 0 && tmp_all_versions[x] <= 3)){
        //check if we've seen this mutation before
        it = std::find(mut.begin(), mut.end(), tmp_all_versions[x]);
        //if we haven't seen it, add it to our repitore
        if(it == mut.end()){
          mut.push_back(tmp_all_versions[x]);
          mut_count.push_back(save_read_counts[x]);
          sub_mut_count.push_back(subtract_read_counts[x]);
          index.push_back({x});
        }else{ //we've seen this before, let's go count it
          int loc_index = it - mut.begin();
          mut_count[loc_index] += save_read_counts[x];
          sub_mut_count[loc_index] += subtract_read_counts[x];
          index[loc_index].push_back(x);
        }
      }
    }
    //we found no mutations that didn't match the reference
    if(mut_count.size() == 0){continue;}
    //go back through each of our saved mutations
    for(uint32_t z= 0; z < mut_count.size(); z++){
      if(mut_count[z] < 10){
        continue;
      }
      double amplicon_freq = mut_count[z] / arc;
      //std::cout << "adjusted_read_count " << adjusted_read_count << " arc " << arc << " mut count " << mut_count[z] <<  " subtract read count " << sub_mut_count[z] << std::endl;

      std::string nuc = decoded_nucs(mut[z], dict_decode);
      double global_freq=0;
      std::vector<allele> pos_alleles = all_positions[final_positions[i]].ad;
      for(allele x : pos_alleles){
        if(x.nuc == nuc){
          global_freq = x.depth / all_positions[final_positions[i]].depth;
        }
      }
      if(abs(global_freq-amplicon_freq) >= 0.20){
        //adjusted_read_count -= mut_count[z];
        suspect_positions.push_back(final_positions[i]);
        //primer_issue = true;
        if(primer_issue){

          double new_total_depth = 0;
          double new_global_depth = 0;
          for(uint32_t y=0; y < pos_alleles.size(); y++){ 
            if(pos_alleles[y].nuc == nuc && pos_alleles[y].ref == false){
              std::cout << "check primer binding " << final_positions[i] << " global_freq " << global_freq << " amplicon_freq " << amplicon_freq << " mut type " << mut[z] << std::endl;
             pos_alleles[y].depth -= sub_mut_count[z];
              if(pos_alleles[y].depth < 0){
                pos_alleles[y].depth = 0;
              }
              new_global_depth = pos_alleles[y].depth;
            }
            new_total_depth += pos_alleles[y].depth;
          }
         double new_global_freq = new_global_depth / new_total_depth;
         if(abs(new_global_freq-amplicon_freq) < 0.20){
           frequencies.push_back(new_global_freq);
         }
         //print_allele_depths(pos_alleles);
         all_positions[final_positions[i]].depth = new_total_depth;
         all_positions[final_positions[i]].ad = pos_alleles;
         masked_positions.push_back(final_positions[i]);
         masked_alleles.push_back(mut[z]);
        }
      } else{
        frequencies.push_back(amplicon_freq);
      }
    }
  }
}

std::vector<double> create_frequency_matrix(IntervalTree &amplicons, std::vector<position> &all_positions, std::vector<primer> primers, std::string output_primer, std::vector<uint32_t> &masked_positions, std::vector<int> &masked_alleles, std::unordered_map<int, std::string> dict_decode, ref_antd reference, std::string region_){
  /*
   * @param amplicons : data strucuture containing read count, haplotype nt, and positions
   * @param all_positions : vector containing depths / alleles for all pos
   * @param masked_primers : vector containing primers with mismatches
   * @return frequencies : a flat vector containing all haplotype frequencies
   *
   * Function calculates the frequency of unique haplotypes on a per amplicon basis. In
   * the process, handles combining haplotypes that originate from reads that don't overlap.
   * Stores unique haplotypes, the associated positions, and frequency to amplicon object.
   *
   * -2 : match ref
   * -1 : soft clip or not covered
   */ 

  std::vector<double> frequencies; //return place
  std::vector<std::vector<uint32_t>> frequency_positions;
  std::vector<std::vector<int>> frequency_haplotypes;

  ITNode *node = amplicons.iterate_nodes();
  int read_count=0; //total reads in amplicon

  //at the amplicon level
  std::vector<std::vector<uint32_t>> positions;
  std::vector<std::vector<uint32_t>> ranges;
  std::vector<std::vector<int>> haplotypes;

  //individual occurences within an amplicon
  std::vector<uint32_t> position;
  std::vector<uint32_t> range;
  std::vector<int> haplotype;
  std::string lower_primer_name;
  
  //loop through all the amplicons
  while(node != NULL){
    lower_primer_name.clear();
    node = amplicons.iterate_nodes(node->right);
    if(node == NULL){
      break;
    }
    read_count = node->read_count;
    if(read_count == 0){
      continue;
    }
    //search for this amplicon as being one with primer binding mutations
    //could use a find function here instead...
    for(uint32_t i = 0; i < primers.size(); i++){
      if(node->data->low == primers[i].get_start()){
        lower_primer_name = primers[i].get_name();
        break;
      }
    }
    double mut_rev_ratio = 0;
    double mut_for_ratio = 0;
    
    if(node->mut_reverse != 0){
      mut_rev_ratio = node->mut_reverse / node->reverse;
    }
    if(node->mut_forward != 0){
      mut_for_ratio = node->mut_forward / node->forward;
    }
    std::string primer_mismatch_percent = std::to_string(mut_for_ratio) + " " + std::to_string(mut_rev_ratio);

    //remove positions & indels where avg. quality is below 20
    remove_low_quality_nts(node, all_positions, dict_decode);

    positions = node->positions;
    haplotypes = node->haplotypes;
    ranges = node->ranges;
    
    //pool all the positions that have been modified in order to create a table
    std::vector<uint32_t> flattened = flatten(positions);
    std::sort(flattened.begin(), flattened.end());
    std::vector<uint32_t>::iterator ip = std::unique(flattened.begin(), flattened.end());
    flattened.resize(std::distance(flattened.begin(), ip));
    std::vector<std::vector<int>> all_unique_haplotypes; //unique occurences of haplotypes on the amplicon
    std::vector<std::vector<int>> all_haplotypes; //extended version of each haplotype covering all positions

    struct position allele_positions;   

    //loop through all the haplotypes in the amplicon and find unqiue ones
    for(uint32_t i=0; i < positions.size(); i++){
      position = positions[i];
      haplotype = haplotypes[i];
      range = ranges[i];
      //initialize vector assuming positions are covered
      std::vector<int> expanded_haplotypes(flattened.size(), -1);

      //fill out this haplotype with -1 for the things that aren't covered or soft clipped
      for(uint32_t i = 0; i < position.size(); i++){
        //calculate allele frequency
        if(haplotype[i] >= 0){
          allele_positions = all_positions[position[i]];
          double freq = allele_positions.ad[haplotype[i]].depth / allele_positions.depth;
          if(freq <= 0.03 || allele_positions.ad[haplotype[i]].depth < 10){
            haplotype[i] = -1;
          }
        }
        std::vector<uint32_t>::iterator it = std::find(flattened.begin(), flattened.end(), position[i]);
        //find this position in our list of all unique pos
        if(it != flattened.end()){
          int index = it - flattened.begin();
          expanded_haplotypes[index] = haplotype[i];
        }
      }
      all_haplotypes.push_back(expanded_haplotypes);
      //have we seen this haplotype before?
      //think this could be cleaned up using sets?
      std::vector<std::vector<int>>::iterator haplo_it = std::find(all_unique_haplotypes.begin(), all_unique_haplotypes.end(),
          expanded_haplotypes);
      //unseen before
      if(haplo_it == all_unique_haplotypes.end()){
        all_unique_haplotypes.push_back(expanded_haplotypes);
      }
    }
    remove_reference_matches(flattened, all_haplotypes);
    if(flattened.size() == 0){
      continue;
    }
        
    std::vector<std::vector<int>> save_haplotypes; //this is where we have our final things
    std::vector<double> save_read_counts; //save the relative count of each haplotype
    std::vector<double> subtract_read_counts;                                          
    double adjusted_read_count = 0;
    
    //the haplotypes by index to be condensed, first one is SC second is matched
    count_haplotype_occurences(all_haplotypes, save_haplotypes, save_read_counts, adjusted_read_count, flattened, subtract_read_counts);

    //save this info to the amplicon
    std::vector<uint32_t> suspect_positions;
    std::vector<uint32_t>::iterator pos_search;

    if(save_haplotypes.size() == 0){continue;}
    bool flagged = false;
    bool primer_issue = true;
    //returns positions that are improperly represented on the amplicon
    check_primer_binding_issues(flattened, all_positions, save_read_counts, adjusted_read_count, save_haplotypes, suspect_positions, masked_positions, masked_alleles, dict_decode, primer_issue, flagged, subtract_read_counts, reference, region_, frequencies);

    //write the suspect positions to a file    
    std::string suspect_position_string = "";
    for(uint32_t i = 0; i < suspect_positions.size(); i++){
      if(i > 0){
        suspect_position_string += "_";
      }
      suspect_position_string += std::to_string(suspect_positions[i]);
    }
    
    std::cout << "sus pos " << suspect_position_string << std::endl;
    std::vector<uint32_t> good_haplotypes;
    std::vector<uint32_t> bad_haplotypes;
    
    //we iterate the haplotypes and whichever contain 'bad' positions they get skipped
    for(uint32_t i = 0; i < save_haplotypes.size(); i++){
      for(uint32_t x = 0; x < save_haplotypes[i].size(); x++){
        pos_search = std::find(masked_positions.begin(), masked_positions.end(), flattened[x]); 
        //check if this haplotype contains a bad pos
        if(pos_search != masked_positions.end()){
          if(save_haplotypes[i][x] == masked_alleles[pos_search - masked_positions.begin()]){
            bad_haplotypes.push_back(i);        
          }
        }else{
         good_haplotypes.push_back(i);
        }
      }
    }

    //get rid of duplicates
    sort(good_haplotypes.begin(), good_haplotypes.end() );
    good_haplotypes.erase(unique(good_haplotypes.begin(), good_haplotypes.end()), good_haplotypes.end());
    sort(bad_haplotypes.begin(), bad_haplotypes.end() );
    bad_haplotypes.erase(unique(bad_haplotypes.begin(), bad_haplotypes.end()), bad_haplotypes.end()); 
    
    //erase the good if they're in the bad
    auto pred = [&bad_haplotypes](const int& key)->bool{
      return std::find(bad_haplotypes.begin(), bad_haplotypes.end(), key) != bad_haplotypes.end();
    };
    good_haplotypes.erase(std::remove_if(good_haplotypes.begin(), good_haplotypes.end(), pred), good_haplotypes.end()); 
    adjusted_read_count = 0;
    for(uint32_t x : good_haplotypes){
      adjusted_read_count += save_read_counts[x];
    }
      
    if(!flagged){
    for(uint32_t i : good_haplotypes){
      node->final_haplotypes.push_back(save_haplotypes[i]);
      //frequencies.push_back(save_read_counts[i] / adjusted_read_count); //record freuqnecies     
      std::cout << save_read_counts[i] / adjusted_read_count << std::endl;                                                               
      //record the haplotype
      frequency_haplotypes.push_back(save_haplotypes[i]);
      //record the positions                                                        
      frequency_positions.push_back(flattened);
      node->frequency.push_back(save_read_counts[i] / adjusted_read_count);
    }
    node->final_positions = flattened;
    //write suspicious positions because of primer issues to text file
    ofstream file;
    if(!output_primer.empty() && suspect_positions.size() > 0){
      file.open(output_primer, ios_base::app);
        file << lower_primer_name << "\t" << suspect_position_string << "\t" << primer_mismatch_percent <<  "\n";
    }
    if(!output_primer.empty()){
      file.close();
    }
    }
  }
  return(frequencies);
}

//acutal consensus call
void call_consensus_from_vector(std::vector<position> all_positions, std::string seq_id, std::string out_file, uint8_t min_qual, double threshold, uint8_t min_depth, char gap, bool min_coverage_flag, double min_insert_threshold){
  /*
   * Function calls consensus on the sequence give a vector of positons where all 
   * allele depths have been calculated.
   */

  std::ofstream fout((out_file+".fa").c_str());
  std::ofstream tmp_qout((out_file+".qual.txt").c_str());
  char *o = new char[out_file.length() + 1];
  strcpy(o, out_file.c_str());
  if(seq_id.empty()) {
    fout << ">Consensus_" << basename(o) << "_threshold_" << threshold << "_quality_" << (uint16_t) min_qual  <<std::endl;
  } else {
    fout << ">" << seq_id <<std::endl;
  }
  delete [] o;
  int mdepth = 0;
  uint32_t prev_pos = 0, pos = 0;
  ret_t t;
  std::string bases;
  std::string qualities;
  std::vector<allele> ad;
  uint32_t bases_zero_depth = 0, bases_min_depth = 0, total_bases = 0;
  for(position p : all_positions){
    ad = p.ad;
    pos = p.pos;
       
    mdepth = p.depth;
    total_bases++;
    if(prev_pos == 0)   // No -/N before alignment starts
      prev_pos = pos;
    if((pos > prev_pos && min_coverage_flag)){
      fout << std::string((pos - prev_pos) - 1, gap);
      tmp_qout << std::string((pos - prev_pos) - 1, '!'); // ! represents 0 quality score.
    }
    if(mdepth >= min_depth){
      t = get_consensus_allele(ad, min_qual, threshold, gap, min_insert_threshold);
      fout << t.nuc;
      tmp_qout << t.q;
    } else{
      bases_min_depth += 1;
      if (mdepth == 0)
  bases_zero_depth += 1;
      if(min_coverage_flag){
  fout << gap;
  tmp_qout << '!';
      }
    }
    ad.clear();
    prev_pos = pos;
  }
  fout << "\n";     // Add new line character after end of sequence
  tmp_qout << "\n";
  tmp_qout.close();
  fout.close();
  std::cout << "Reference length: " << total_bases << std::endl;
  std::cout << "Positions with 0 depth: " << bases_zero_depth << std::endl;
  std::cout << "Positions with depth below " <<(unsigned) min_depth << ": " << bases_min_depth << std::endl;
}

//entry point for threshold determination
int determine_threshold(std::string bam, std::string ref, std::string bed, std::string pair_info, int32_t primer_offset, double min_insert_threshold, uint8_t min_qual, char gap, double min_depth, bool min_coverage_flag, std::string prefix){
  /*
   * @param bam : path to the bam file
   * @param ref : path to the reference file
   * @param bed : path to the bed file
   * @param pair_info : path to the primer pair .tsv file
   * @param primer_offset : 
   */
  std::unordered_map<std::string, int> dict_encode;
  std::unordered_map<int, std::string> dict_decode;

  dict_encode["A"] = 0;
  dict_encode["S"] = -1;
  dict_encode["C"] = 1;
  dict_encode["G"] = 2;
  dict_encode["T"] = 3;
  dict_encode["*"] = 4;

  dict_decode[0] = "A";
  dict_decode[1] = "C";
  dict_decode[2] = "G";
  dict_decode[3] = "T";
  dict_decode[-1] = "S";
  dict_decode[4] = "*";

  //TODO: pass this as a param 
  bool mask_primer_muts = false;

  //boiler plate
  //generate the .fa consensus file header
  std::string suffix = ".bam";
  std::string seq_id = "Consensus_" + bam.substr(0, bam.length() - suffix.length());
  std::vector<position> all_positions;
  std::string output_amplicon = prefix + ".txt";
  std::vector<primer> primers;
  IntervalTree amplicons;
  int read_counter = 0;
  std::string output_primer = prefix + "_primer_mismatches.txt";
  ofstream file;
  alglib::real_2d_array xy;
  int best_cluster_index = 0;
  double best_sil_score = 0;
  double threshold = 0;  
  uint32_t max_n = 6; //max number of clusters to attempt
  std::vector<cluster> all_cluster_results;
  int tmp_cluster_index;
  double tmp_thresh;
  double sum_centers = 0;
  int i = 0;
  std::string cluster_filename = prefix + "_cluster_results.txt";

  //preset the alleles to save time later
  std::vector<std::string> basic_nts = {"A", "C", "G", "T"};
  std::vector<allele> basic_alleles;
  for(std::string nt : basic_nts){
    allele new_allele;
    new_allele.nuc = nt;
    new_allele.depth = 0;
    new_allele.tmp_mean_qual=0;
    basic_alleles.push_back(new_allele);
  }

  //populate primer, and primer pairs
  primers = populate_from_file(bed, primer_offset);
  amplicons = populate_amplicons(pair_info, primers);
  
  samFile *in = hts_open(bam.c_str(), "r");  
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  bam_hdr_t *header = sam_hdr_read(in);

  uint32_t *ref_length = header->target_len;
  bam1_t *aln = bam_init1();
  hts_itr_t *iter = NULL;
  std::string region_;
  //region refers to reference
  region_.assign(header->target_name[0]);
  iter = sam_itr_querys(idx, header, region_.c_str());
  
  //load the reference
  ref_antd reference(ref);
  reference.set_seq(region_);

  std::string ref_nuc;
  //fill md tag TODO
  //bam_fillmd1_core(const char *ref_name, aln, char *ref, int flag, int max_nm)
  for(uint32_t i=0; i <= *ref_length; i++){
    position new_position;
    new_position.pos = i;
    new_position.depth = 0;
    new_position.ad = basic_alleles;    
    ref_nuc += reference.get_base(i, region_);
    for(uint32_t i = 0; i < new_position.ad.size(); i++){
      if(new_position.ad[i].nuc == ref_nuc){
        new_position.ad[i].ref = true;
        break;
      }
    }
    ref_nuc.clear();
    all_positions.push_back(new_position);
  }

  //this iterates over the reads and assigns them to an amplicon
  while(sam_itr_next(in, iter, aln) >= 0) {
    //TODO: reimagine this as a percent
    if(read_counter % 100000 == 0){
      std::cout << read_counter << " reads processed." << std::endl;
    }
    if(read_counter < 1400000 || read_counter > 1500000){
      read_counter += 1;
      continue;
    }
    read_counter += 1;
    iterate_reads(aln, amplicons, all_positions, reference, region_, dict_encode, dict_decode);
  }
  std::cout << "end read processing." << std::endl; 
 
  //calculate the mean quality, might not print properly but saves properly
  for(uint32_t i = 0; i < all_positions.size(); i++){
    for(uint32_t x = 0; x < all_positions[i].ad.size(); x++){
      allele al = all_positions[i].ad[x];
      uint8_t tmp_qual = (uint8_t) floor((al.tmp_mean_qual / al.depth)+0.5);
      all_positions[i].ad[x].mean_qual = tmp_qual;
    }  
  }

  file.open(output_primer, ios_base::app);
  file << "lower_primer_name" << "\t" << "suspect_positions" << "\t" << "mutation_percent_primer" << "\n";
  file.close(); 
  std::vector<uint32_t> masked_positions;
  std::vector<int> masked_alleles;
  
  print_allele_depths(all_positions[23604].ad);
  //extract those reads into a format useable in the clustering
  std::vector<double> all_frequencies = create_frequency_matrix(amplicons, all_positions, primers, output_primer, masked_positions, masked_alleles, dict_decode, reference, region_);
  reference.remove_seq(); 
  print_allele_depths(all_positions[23604].ad);
  if(all_frequencies.size() < 2){
    return(0);
  }
  //remove perfect 1 haplotypes
  all_frequencies.erase(std::remove_if(
    all_frequencies.begin(), all_frequencies.end(),
    [](double& x) { 
        return(x==1 || x < 0.001);
    }), all_frequencies.end());

  //if we're going to masked the mutations, remove them here
  if(mask_primer_muts){
    for(uint32_t i=0; i < masked_positions.size(); i++){
      position tmp_pos = all_positions[masked_positions[i]];
      std::string tmp_nuc = decoded_nucs(masked_alleles[i], dict_decode);
      for(uint32_t x = 0; x < tmp_pos.ad.size(); x++){
        if(tmp_pos.ad[x].nuc == tmp_nuc){
          std::cout << masked_positions[i] << " " << masked_alleles[i] << " " << tmp_nuc << std::endl;
          tmp_pos.depth -= tmp_pos.ad[x].depth;
          tmp_pos.ad[x].depth = 0;
          all_positions[masked_positions[i]] = tmp_pos;
        }     
      }
    }
  }
  file.open(output_amplicon, ios_base::app);
  file << "lower_primer\tupper_primer\tread_count\tpositions\tfrequencies\thaplotypes\tnumber_haplotypes\n";
  file.close();
  amplicons.dump_amplicon_summary(output_amplicon);
  //reshape it into a real 2d array for alglib
  xy.setlength(all_frequencies.size(), 1);
  for(uint32_t i=0; i < all_frequencies.size(); i++){
    xy(i,0) = all_frequencies[i];
  }

  std::sort(all_frequencies.begin(), all_frequencies.end());
  auto uniq = std::unique(all_frequencies.begin(), all_frequencies.end()) - all_frequencies.begin();
  //if we have fewer than 6 points, we can only have that many clusters
  if(uniq < max_n){
    max_n = uniq;
  }
  //call kmeans clustering
  for (uint32_t n = 2; n < max_n; n++){
    //call kmeans clustering
    cluster cluster_results; //reset the results
    k_means(n, xy, cluster_results);
    all_cluster_results.push_back(cluster_results);
  }
  //open the file to save clustering results
  file.open(cluster_filename, ios_base::app);
  file << "sil_score\tadjusted_sil_score\tcluster_centers\tn_clusters\tthreshold\n";
  //dump additional information to a file such as (1) cluster values (2) cluster centers
  for(cluster x : all_cluster_results){
    sum_centers = 0;
    //threshold for if we chose this cluster
    tmp_cluster_index = std::max_element(x.centers.begin(), x.centers.end()) - x.centers.begin();
    tmp_thresh = x.cluster_bounds[tmp_cluster_index][0] - 0.01;
    file << x.sil_score << "\t";
    for(double cent : x.centers){
      sum_centers += cent;
    }
    file << x.sil_score - (abs(sum_centers -1)) << "\t";
    if((x.sil_score - abs(sum_centers -1))  > best_sil_score){
      best_sil_score = x.sil_score - abs(sum_centers-1);
      best_cluster_index = i;
   }
   for(double c : x.centers){
      file << c << "_";
    }
    file << "\t";
    file << x.n_clusters << "\t";
    file << tmp_thresh << "\n";
    i++;
  }
  file.close();
  cluster choice_cluster = all_cluster_results[best_cluster_index];
  //find the largest cluster center
  int largest_cluster_index = std::max_element(choice_cluster.centers.begin(), choice_cluster.centers.end()) - choice_cluster.centers.begin();
  //this marks the lower bound of the largest cluster
  threshold = choice_cluster.cluster_bounds[largest_cluster_index][0] - 0.01;

  //determine whether or not we should call consensus
  /*if(choice_cluster.sil_score <= 0.80 || threshold <= 0.5){
    return(0);
  }*/
   
  //call consensus
  call_consensus_from_vector(all_positions, seq_id, prefix, min_qual, threshold, min_depth, gap, min_coverage_flag, min_insert_threshold);
  
  file.open("frequencies_corrected.txt", ios_base::app);
  file << "position\tfrequency\tnuc\ttotal_depth\n";
 //dump out frequency information for my own personal use
  for(position p : all_positions){
    for(allele a : p.ad){
      if(a.depth > 0){
        file << p.pos << "\t" << a.depth/p.depth << "\t" << a.nuc << "\t"<< p.depth <<"\n"; 
      }
    }
  }
  file.close();

  return 0;
}

