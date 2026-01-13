#include <vector>
#include <fstream>
#ifndef gmm
#define gmm

struct gaussian_mixture_model {
  std::vector<std::vector<double>> prob_matrix;
  uint32_t n;
  uint32_t lower_n;
  double var_floor;
  double bic;
  double log_likelihood;
  std::vector<double> means;
  std::vector<double> hefts;
  std::vector<double> dcovs;
  std::vector<std::vector<double>> clusters;
  std::vector<double> cluster_std_devs;
  std::vector<double> cluster_probabilities;
};

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  uint32_t total_depth;
  uint32_t gapped_depth;
  double qual;
  double freq;
  double gapped_freq;
  bool version_1_var=false;
  double std_dev;

  //number corresponding the the amplicons covering this position
  std::vector<uint32_t> amplicon_numbers;
  //frequencies of this variants on each amplicon
  std::vector<double> freq_numbers;
  //per amplicon frequency assignments to clusters
  std::vector<uint32_t> freq_assignments;
  //the consensus sequence this variant is assigned to
  std::vector<uint32_t> consensus_numbers;
  //if this cluster is fully resolveable or not
  bool resolved=true;

  //for these true means flagged as problematic
  bool vague_component_assignment=false; //variant cannot be distinguished between two components
  bool amplicon_flux=false; //fluctuation frequency across amplicons
  bool amplicon_masked=false; //masked due to another variant experiencing flux
  bool depth_flag=false; //depth is below the threshold
  bool qual_flag=false; //quality is below threshold
  bool outside_freq_range=false; //outside of useful frequency range for model
  bool include_clustering=true; //here we flag the later positions of deletions

  std::vector<double> marginal_posterior_probabilities;
  int assigned_component = -1;
};

void perm_generator(int n, int k, std::vector<std::vector<uint32_t>> &possible_permutations);
void split(std::string &s, char delim, std::vector<std::string> &elems);
int gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, std::string ref, double default_threshold);
void parse_internal_variants(std::string filename, std::vector<variant> &base_variants, uint32_t depth_cutoff, uint32_t round_val, uint8_t quality_threshold);
uint32_t smallest_value_index(std::vector<double> values);
std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>>& input_vector);
double calculate_mean(const std::vector<double>& data);
void set_freq_range_flags(std::vector<variant> &variants, double lower_bound, double upper_bound, bool advanced);
void set_deletion_flags(std::vector<variant> &variants, double lower_bound);
void noise_resampler(uint32_t n, uint32_t index, std::vector<std::vector<uint32_t>> &possible_permutations, uint32_t amount_resample);
#endif