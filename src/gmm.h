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
  bool half_normal_upper = false;
  bool half_normal_lower = false;
  bool position_half_normal_upper = false;

  //identifiers of the amplicons covering this position, "low-high"
  std::vector<std::string> amplicon_ids;
  //frequencies of this variants on each amplicon
  std::vector<double> freq_numbers;
  //depth of each amplicon covering this position
  std::vector<uint32_t> amplicon_depths;
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
  bool position_masked=false; //this position fluctuates in frequency across amplicons
  bool amplicon_masked=false; //an amplicon covering this position fluctuates elsewhere
  bool depth_flag=false; //depth is below the threshold
  bool qual_flag=false; //quality is below threshold
  bool outside_freq_range=false; //outside of useful frequency range for model
  bool cluster_outlier=false; //is an outlier for the cluster assigned
  bool wide_sd=false; //further than WIDE_SD_THRESHOLD sds from the cluster it was assigned to
  bool overlapped_deletion=false; //minor deletion overlapping a more abundant deletion at the same site
  bool imbalance=false;
  bool position_conflict=false; //multiple variants at this position assigned to the same cluster
  std::vector<double> probabilities;

};

extern const double DEFAULT_AMPLICON_STDEV;

//a variant further than this many sds from its assigned cluster is not a credible
//member of it, so the position is called N rather than given a base
static constexpr double WIDE_SD_THRESHOLD = 3.0;

void flag_wide_sd_variants(std::vector<variant> &variants, const std::vector<double> &eff_means,
                           const std::vector<double> &eff_vars, const std::vector<double> &unrefined_means);

void split(std::string &s, char delim, std::vector<std::string> &elems);
std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, std::vector<double> &solution, std::vector<double> &means, double default_threshold, uint32_t n, double invariant_threshold, double covariance_prior, double mean_precision_prior, double half_normal_covariance_prior, double min_cluster_fraction = 0.10, uint32_t min_cluster_points = 0, double amplicon_stdev = DEFAULT_AMPLICON_STDEV);
void parse_internal_variants(std::string filename, std::vector<variant> &base_variants, uint32_t depth_cutoff, uint32_t round_val, uint8_t quality_threshold, double invariant_threshold);
std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>>& input_vector);
void set_freq_range_flags(std::vector<variant> &variants, double lower_bound, double upper_bound, bool advanced);
void set_deletion_flags(std::vector<variant> &variants, double lower_bound, double invariant_lower_bound);
void flag_amplicon_variation(std::vector<variant> &variants, double amplicon_stdev = DEFAULT_AMPLICON_STDEV);
void propagate_amplicon_masking(std::vector<variant> &variants);
void rewrite_position_masking(std::vector<variant> &variants);
#endif