#include <htslib/faidx.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "allele_functions.h"
#include "parse_gff.h"

#ifndef ref_seq
#define ref_seq

const char UNKNOWN_BASE = 'N';
extern const unsigned char comp_base[256];

class ref_antd {
 public:
  ref_antd(std::string ref_path);
  ref_antd(std::string ref_path, std::string gff_path);
  ~ref_antd();
  char get_base(int64_t pos, std::string region);
  int64_t get_length(std::string region);
  void complement_codon(char* codon);
  void reverse_complement_codon(char* codon);
  int add_gff(std::string path);
  int set_index(std::string path);

  void set_seq(std::string &region);

  int codon_aa_stream(std::string region, std::ostringstream &line_stream,
                      std::ofstream &fout, int64_t pos, char alt);
  char *get_codon(int64_t pos, std::string region, gff3_feature feature);
  char *get_codon(int64_t pos, std::string region, gff3_feature feature,
                  char alt);
  char *get_codon(int64_t pos, std::string region, const cds_group &group);
  std::vector<cds_group> query_cds_groups(int64_t pos);
  const std::vector<cds_group> &get_cds_groups() const;
  std::vector<gff3_feature> get_gff_features();

 private:
  char *seq;
  gff3 gff;
  faidx_t *fai;
  std::string region;
};

#endif
