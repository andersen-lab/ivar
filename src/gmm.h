#include <vector>
#include <fstream>
#include "./include/armadillo"
#ifndef gmm
#define gmm

struct kmeans_model {
  std::vector<std::vector<double>> clusters; //stored assigned clusters
  uint32_t n; //number of clusters
  std::vector<double> means; //centroids
};


struct gaussian_mixture_model {
  std::vector<std::vector<double>> prob_matrix;
  uint32_t n;
  double var_floor;
  std::vector<double> means;
  std::vector<double> hefts; 
  std::vector<double> dcovs;
  double aic;
  arma::gmm_diag model;
  std::vector<std::vector<double>> clusters;
  std::vector<double> cluster_std_devs;
};

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  double qual;
  double freq;
  double gapped_freq = 0;
  int cluster_assigned = -1;
  bool version_1_var=false;  
  double std_dev;
  //number corresponding the the amplicons covering this position
  std::vector<uint32_t> amplicon_numbers;
  //frequencies of this variants on each amplicon
  std::vector<double> freq_numbers;
  //the consensus sequence this variant is assigned to
  std::vector<uint32_t> consensus_numbers;
  //if this cluster is fully resolveable or not
  bool resolved=true;

  //for these true means flagged as problematic
  bool vague_assignment=false; //cannot be distinguished between two groups
  bool amplicon_flux=false; //fluctuation frequency across amplicons
  bool amplicon_masked=false; //masked due to another variant experiencing flu
  bool primer_masked=false; //mutation in primer binding region of overlapped amplicon
  bool depth_flag=false; //depth is below the threshold                  
  bool low_prob_flag=false;
  bool qual_flag=false;
  bool outside_freq_range=false; //outside of useful frequency range for model                 
  bool cluster_outlier=false; //is an outlier for the cluster assigned
  std::vector<double> probabilities;

};
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<variant> gmm_model(std::string prefix, std::string output_prefix, uint32_t min_depth, uint8_t min_qual, std::vector<double> &solution, std::vector<double> &means);
void parse_internal_variants(std::string filename, std::vector<variant> &base_variants, uint32_t depth_cutoff, uint32_t round_val, uint8_t quality_threshold);
uint32_t smallest_value_index(std::vector<double> values);
std::vector<std::vector<double>> transpose_vector(const std::vector<std::vector<double>>& input_vector);
void assign_variants_simple(std::vector<variant> &variants, std::vector<std::vector<double>> prob_matrix, uint32_t index, uint32_t lower_n);
gaussian_mixture_model retrain_model(uint32_t n, arma::mat data, std::vector<variant> &variants, uint32_t lower_n, double var_floor);
void assign_clusters(std::vector<variant> &variants, gaussian_mixture_model gmodel, uint32_t lower_n);
double calculate_mean(const std::vector<double>& data);
void calculate_reference_frequency(std::vector<variant> &variants, std::string filename, uint32_t depth_cutoff, double lower_bound, double upper_bound);
kmeans_model train_model(uint32_t n, arma::mat data, bool error);
void set_freq_range_flags(std::vector<variant> &variants, double lower_bound, double upper_bound);
#endif
