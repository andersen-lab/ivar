#include <vector>
#include <fstream>
#ifndef gmm
#define gmm

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  uint32_t total_depth;
  uint32_t gapped_depth;
  double qual;
  double freq;
  double gapped_freq = 0;
  double logit = 0;
  int cluster_assigned = -1;
  bool version_1_var=false;
  double std_dev;
  bool half_normal_upper = false;
  bool half_normal_lower = false;
  bool position_half_normal_upper = false;

  //number corresponding the the amplicons covering this position
  std::vector<uint32_t> amplicon_numbers;
  //frequencies of this variants on each amplicon
  std::vector<double> freq_numbers;
  //per amplicon frequency assignments to clusters
  std::vector<uint32_t> freq_assignments;
  //the consensus sequence this variant is assigned to
  std::vector<uint32_t> consensus_numbers;
  //consensus genomes for which this variant's peak has multiple explanations
  std::vector<uint32_t> ambiguous_numbers;
  //if this cluster is fully resolveable or not
  bool resolved=true;

  bool assigned_deletion=false;

  //for these true means flagged as problematic
  bool amplicon_flux=false; //fluctuation frequency across amplicons
  bool amplicon_masked=false; //masked due to another variant experiencing flux
  bool primer_masked=false; //mutation in primer binding region of overlapped amplicon
  bool depth_flag=false; //depth is below the threshold
  bool qual_flag=false; //quality is below threshold
  bool outside_freq_range=false; //outside of useful frequency range for model
  bool cluster_outlier=false; //is an outlier for the cluster assigned
  bool overlapped_deletion=false; //minor deletion overlapping a more abundant deletion at the same site
  bool imbalance=false;
  bool position_conflict=false; //multiple variants at this position assigned to the same cluster
  std::vector<double> probabilities;

};

void split(std::string &s, char delim, std::vector<std::string> &elems);
std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, std::vector<double> &solution, std::vector<double> &means, double default_threshold, uint32_t n, double invariant_threshold, double covariance_prior, double mean_precision_prior, double min_cluster_fraction = 0.10);
void parse_internal_variants(std::string filename, std::vector<variant> &base_variants, uint32_t depth_cutoff, uint32_t round_val, uint8_t quality_threshold, double invariant_threshold);
std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>>& input_vector);
void set_freq_range_flags(std::vector<variant> &variants, double lower_bound, double upper_bound, bool advanced);
void set_deletion_flags(std::vector<variant> &variants, double lower_bound, double invariant_lower_bound);
void reset_variants_info(std::vector<variant> &variants);
void rewrite_position_masking(std::vector<variant> &variants);
#endif