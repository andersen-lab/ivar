#include <iostream>
#include <sstream>
#include <fstream>

#include <htslib/faidx.h>
#include "parse_gff.h"
#include "allele_functions.h"

#ifndef ref_seq
#define ref_seq

const char UNKNOWN_BASE = 'N';

class ref_antd{
public:
  ref_antd(std::string ref_path);
  ref_antd(std::string ref_path, std::string gff_path);
  ~ref_antd();
  char get_base(int64_t pos, std::string region);
  void set_seq(std::string region);
  void remove_seq();
  int add_gff(std::string path);
  int add_seq(std::string path);
  int codon_aa_stream(std::string region, std::ostringstream &line_stream, std::ofstream &fout, int64_t pos, char alt);
  char* get_codon(int64_t pos, std::string region, gff3_feature feature);
  char* get_codon(int64_t pos, std::string region, gff3_feature feature, char alt);
  std::vector<gff3_feature> get_gff_features();

private:
  char *seq;
  gff3 gff;
  faidx_t *fai;
  std::string region;
};

#endif
