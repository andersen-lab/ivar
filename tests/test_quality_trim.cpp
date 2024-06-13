#include <iostream>

#include "../src/trim_primer_quality.h"
#include "htslib/sam.h"

int main() {
  int success = 0;
  std::string bam = "../data/test.sorted.bam";
  std::string region_;
  samFile *in = hts_open(bam.c_str(), "r");
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  bam_hdr_t *header = sam_hdr_read(in);
  region_.assign(header->target_name[0]);
  std::string temp(header->text);
  hts_itr_t *iter = NULL;
  iter = sam_itr_querys(idx, header, region_.c_str());
  bam1_t *aln = bam_init1();
  cigar_ t;
  int lengths[6] = {150, 99, 99, 24, 146, 144}, ctr = 0;
  int start_pos_rev[6] = {19, 114, 209, 325, 199, 231};
  while (sam_itr_next(in, iter, aln) >= 0) {
    t = quality_trim(aln, 20, 4);
    std::cout << bam_get_qname(aln) << std::endl;
    std::cout << "POS: " << t.start_pos
              << " Length: " << bam_cigar2rlen(t.nlength, t.cigar) << std::endl;
    if (t.start_pos != start_pos_rev[ctr]) {
      success = -1;
    }
    if (bam_cigar2rlen(t.nlength, t.cigar) != lengths[ctr]) {
      success = -1;
    }
    free_cigar(t);
    ctr++;
  }

  bam_destroy1(aln);
  bam_itr_destroy(iter);
  sam_hdr_destroy(header);
  hts_idx_destroy(idx);
  hts_close(in);
  return success;
}
