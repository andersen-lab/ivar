#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

class consensus_sequence {
  std::string seq_name;
  uint32_t seq_length;
  std::string sequence;
  std::vector<std::vector<variant>> variant_records;

  public:
    consensus_sequence(uint32_t max_position) {
      seq_length = max_position;
      sequence = std::string(max_position, 'N');
      variant_records = std::vector<std::vector<variant>>(max_position);
    }

    void add_variant(uint32_t position, const variant &v) {
      variant_records[position-1].push_back(v);
    }
};

std::string trim_leading_ambiguities(std::string sequence, uint32_t min_position);
void cluster_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means);
#endif
