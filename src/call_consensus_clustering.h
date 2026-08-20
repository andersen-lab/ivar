#include "gmm.h"
#ifndef call_consensus_clustering
#define call_consensus_clustering

class consensus_sequence {
  std::string seq_name;
  uint32_t seq_length;
  std::vector<std::string> sequence;
  std::vector<std::vector<variant>> variant_records;
  std::vector<bool> position_ambiguous;

  public:
    consensus_sequence(uint32_t max_position) {
      seq_length = max_position;
      sequence = std::vector<std::string>(max_position, "N");
      variant_records = std::vector<std::vector<variant>>(max_position);
      position_ambiguous = std::vector<bool>(max_position, false);
    }

    void add_variant(uint32_t position, const variant &v) {
      variant_records[position-1].push_back(v);
    }

    void mark_ambiguous(uint32_t position) {
      position_ambiguous[position-1] = true;
    }

    void set_seq_name(std::string name) {
      seq_name = name;
    }

    std::string get_base(uint32_t position) const {
      return sequence[position-1];
    }

    void get_consensus(uint32_t n);
    void get_majority_consensus(double threshold);
    void process_variant_assignments();
    void write_consensus_to_file(std::string consensus_filename);
};

void cluster_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means);
void call_majority_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold);
void assign_variants_position(std::vector<variant> &variants, std::vector<consensus_sequence> &all_consensus_seqs);
#endif
