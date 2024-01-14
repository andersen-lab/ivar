#ifndef gmm
#define gmm

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  float qual;
  float freq;
  //for these true means flagged as problematic
  bool vague_assignment; //cannot be distinguished between two groups
  bool amplicon_flux; //fluctuatin frequency across amplicons
  bool depth_flag; //depth is below the threshold                  
  bool outside_freq_range; //outside of useful frequency range for model                 
  std::vector<double> probabilities;

};

int gmm_model(std::string prefix);
#endif
