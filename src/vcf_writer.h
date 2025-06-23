#include <fstream>
#include <iostream>
#include <vector>

#include "allele_functions.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "ref_seq.h"

class vcf_writer {
  vcfFile *file;
  ref_antd *ref;
  bcf_hdr_t *hdr;
  std::string region;
  std::string sample_name;
  int write_record(uint32_t pos, std::vector<allele> &alleles, char ref);
  int init_header();

 public:
  int init(char _mode, std::string fname, std::string region,
           std::string sample_name, std::string ref_path);
};
