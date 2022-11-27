#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "./alglib/ap.h"
#include "./alglib/dataanalysis.h"
#include "interval_tree.h"
#include "./alglib/stdafx.h"
#include "primer_bed.h"

#ifndef clustering
#define clustering
class read_linkage {
  public:
    std::string read_name;
    std::vector<uint32_t> positions;
    std::vector<int> haplotype;

};
class cluster {
  public:
    int n_clusters;
    double sil_score; //average of per sample sil score
    std::vector<double> sil_scores; //per sample                                       
    std::vector<double> centers;
    std::vector<std::vector<double>> sorted_points; //all points sorted by center they belong to
    std::vector<std::vector<double>> cluster_bounds; //the upper and lower bounds of each cluster
};

//zipping and unzipping functions for haplotypes && position
void zip(const std::vector<int> &haplotypes, const std::vector<uint32_t> &positions,std::vector<std::pair<uint32_t,int>> &zipped);
void unzip(const std::vector<std::pair<uint32_t, int>> &zipped, std::vector<int> &haplotypes, std::vector<uint32_t> &positions);
void print_cluster_info(cluster cluster_results);
void k_means(int n_clusters, alglib::real_2d_array xy, cluster &cluster_results);
std::vector<double> calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters);
double cluster_point_distances(alglib::real_2d_array X, alglib::kmeansreport rep, double point, double center, int n_clusters);
void calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep,int n_clusters, cluster &cluster_results);
std::vector<double> create_frequency_matrix(IntervalTree &amplicons, std::vector<position> all_positions, std::vector<primer> masked_primers, std::string output_primer);
int determine_threshold(std::string bam, std::string ref, std::string bed, std::string pair_info, int32_t primer_offset, double min_insert_threshold, uint8_t min_qual, char gap, double min_depth, bool min_coverage_flag, std::string prefix);
void parse_md_tag(uint8_t *aux, std::vector<int> &haplotypes, std::vector<uint32_t> &positions,
    uint32_t abs_start_pos, std::vector<position> &all_positions,
    uint8_t *seq, uint32_t length, uint32_t correction_factor, uint32_t abs_end_pos,
    std::vector<uint32_t> ignore_positions, bool reverse);
std::string decoded_nucs(int tmp, std::unordered_map<int, std::string> dict_decode);
std::vector<std::vector<int>> transpose(const std::vector<std::vector<int>> data);
#endif
