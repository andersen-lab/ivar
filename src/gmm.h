#ifndef gmm
#define gmm

struct variant {
  uint32_t position;
  std::string nuc;
  uint32_t depth;
  float qual;
  float freq;
  bool vague_assignment;
}

int gmm_model(std::string prefix);
#endif
