#include <vector>
#include <fstream>
#include "./include/armadillo"
#ifndef gmm
#define gmm

struct gaussian_mixture_model {
  std::vector<std::vector<double>> prob_matrix;
  uint32_t n;
  double var_floor;
  std::vector<double> means;
  std::vector<double> hefts; 
  std::vector<double> dcovs;
  arma::gmm_diag model;
};

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  float qual;
  float freq;
  float gapped_freq = 0;
  float transformed_freq;
  float transformed_gap_freq;
  int cluster_assigned = -1;
  bool version_1_var=false;  
  //for these true means flagged as problematic
  bool vague_assignment=false; //cannot be distinguished between two groups
  bool amplicon_flux=false; //fluctuation frequency across amplicons
  bool amplicon_masked=false; //masked due to another variant experiencing flu
  bool primer_masked=false; //mutation in primer binding region of overlapped amplicon
  bool depth_flag=false; //depth is below the threshold                  
  bool low_prob_flag=false;
  bool del_flag=false;
  bool pos_del_flag=false;
  bool qual_flag=false;
  bool outside_freq_range=false; //outside of useful frequency range for model                 
  bool gap_outside_freq_range=false;
  bool cluster_outlier=false; //is an outlier for the cluster assigned
  std::vector<double> probabilities;

};
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<variant> gmm_model(std::string prefix, std::string output_prefix);
void parse_internal_variants(std::string filename, std::vector<variant> &base_variants, uint32_t depth_cutoff, float lower_bound, float upper_bound, std::vector<uint32_t> deletion_positions, uint32_t round_val);
std::vector<uint32_t> find_deletion_positions(std::string filename, uint32_t depth_cutoff, float lower_bound, float upper_bound, uint32_t round_val);
std::vector<uint32_t> find_low_quality_positions(std::string filename, uint32_t depth_cutoff, float lower_bound, float upper_bound, float quality_threshold, uint32_t round_val);
std::vector<std::vector<double>> solve_possible_solutions(std::vector<float> tmp_means, double error);
uint32_t smallest_value_index(std::vector<double> values);
std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>>& input_vector);
void assign_variants_simple(std::vector<variant> &variants, std::vector<std::vector<double>> prob_matrix, uint32_t index, std::vector<double> means, uint32_t lower_n);
gaussian_mixture_model retrain_model(uint32_t n, arma::mat data, std::vector<variant> variants, uint32_t lower_n, double var_floor);
void assign_clusters(std::vector<variant> &variants, gaussian_mixture_model gmodel, uint32_t lower_n);
double calculate_mean(const std::vector<double>& data);
void calculate_reference_frequency(std::vector<variant> &variants, std::string filename, uint32_t depth_cutoff, float lower_bound, float upper_bound, std::vector<uint32_t> deletion_positions);
#endif
