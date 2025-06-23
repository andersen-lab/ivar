#include <fstream>
#include <iostream>
#include <vector>

#include "../src/remove_reads_from_amplicon.h"
#include "htslib/sam.h"

int main() {
  int num_success = 0, num_tests = 1;
  std::vector<std::string> amp;
  std::ifstream fin("../data/test.masked_primer_indices.txt");
  std::string s, region = "Consensus_ZI-27_threshold_0_quality_20";
  while (getline(fin, s, '\t')) {
    amp.push_back(s);
  }
  rmv_reads_from_amplicon(
      "../data/test.trimmed.sorted.bam", "", "../data/test.trimmed.masked", amp,
      "../data/test.bed",
      "@PG\tID:ivar-removereads\tPN:ivar\tVN:1.0.0\tCL:ivar removereads -i "
      "../data/test.trimmed.sorted.bam -p ../data/test.trimmed.masked -t "
      "../data/test.masked_primer_indices.txt -b ../data/test.bed\n\0");
  std::string out_file = "../data/test.trimmed.masked.bam";
  samFile *in = hts_open(out_file.c_str(), "r");
  bam_hdr_t *header = sam_hdr_read(in);
  hts_itr_t *iter = NULL;
  if (sam_index_build2(out_file.c_str(), 0, 0) < 0) {
    std::cerr << "Failed to build index" << std::endl;
    return 1;
  }
  hts_idx_t *idx = sam_index_load(in, out_file.c_str());
  iter = sam_itr_querys(idx, header, region.c_str());
  bam1_t *aln = bam_init1();
  bool w;
  while (sam_itr_next(in, iter, aln) >= 0) {
    uint8_t *a = bam_aux_get(aln, "XA");
    w = true;
    if (a != 0) {
      if (bam_aux2i(a) == 2 || bam_aux2i(a) == 4 || bam_aux2i(a) == 6 ||
          bam_aux2i(a) == 1 || bam_aux2i(a) == 3 || bam_aux2i(a) == 5) {
        w = false;
      }
    }
    if (!w) break;
  }
  if (w) num_success += 1;
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return (num_success == num_tests) ? 0 : 1;
}
