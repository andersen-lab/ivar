#include <vector>
#ifndef gmm
#define gmm

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  float qual;
  float freq;
  bool is_ref=false; 
  int cluster_assigned = -1;

  //for these true means flagged as problematic
  bool vague_assignment=false; //cannot be distinguished between two groups
  bool amplicon_flux=false; //fluctuation frequency across amplicons
  bool depth_flag=false; //depth is below the threshold                  
  bool outside_freq_range=false; //outside of useful frequency range for model                 
  bool cluster_outlier=false; //is an outlier for the cluster assigned
  std::vector<double> probabilities;

};

int gmm_model(std::string prefix);
void parse_internal_variants(std::string filename, std::vector<variant> &variants, uint32_t depth_cutoff, float lower_bound, float upper_bound);
#endif
