#include <iostream>
#include <vector>

#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/trim_primer_quality.h"
#include "htslib/sam.h"

int main() {
  int num_tests = 1;
  int success = 0;
  std::string bam =
      "../data/primer/pair_test/test_primer_pair_trim.untrimmed.bam";
  std::string pair_info = "../data/primer_pair_test/test_primer_pair_trim.tsv"

      // populate primers and amplicons from pair file
      std::vector<primer>
          primers = populate_from_file(
              "../data/primer_pair_test/test_primer_pair.bed");
  int max_primer_len = 0;
  max_primer_len = get_bigger_primer(primers);
  std::string region_;
  // open bam file and build and index if needed
  samFile *in = hts_open(bam.c_str(), "r");
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  if (idx == NULL) {
    if (sam_index_build2(bam.c_str(), 0, 0) < 0) {
      std::cerr << ("Unable to open BAM/SAM index.") << std::endl;
      return -1;
    } else {
      idx = sam_index_load(in, bam.c_str());
      if (idx == NULL) {
        std::cerr << "Unable to create BAM/SAM index." << std::endl;
        return -1;
      }
    }
  }
  // read in header
  bam_hdr_t *header = sam_hdr_read(in);
  region_.assign(header->target_name[0]);
  std::string temp(header->text);
  hts_itr_t *iter = NULL;
  iter = sam_itr_querys(idx, header, region_.c_str());
  bam1_t *aln = bam_init1();
  cigar_ t;
  uint32_t *cigar;
  int primer_ctr = 0;
  // all this needs to be modified
  int primer_indices[] = {0, 0, 7, 7, 6};
  uint8_t cigar_flag[5][11] = {
      {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CMATCH},
      {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP,
       BAM_CSOFT_CLIP, BAM_CDEL, BAM_CPAD, BAM_CDEL, BAM_CMATCH},
      {BAM_CMATCH, BAM_CPAD, BAM_CDEL, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP,
       BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP},
      {BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP,
       BAM_CSOFT_CLIP, BAM_CSOFT_CLIP},
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CDEL, BAM_CMATCH, BAM_CPAD, BAM_CMATCH,
       BAM_CPAD, BAM_CMATCH}};
  uint32_t read_start_pos[5] = {
      30, 29, 231, 249, 274  // also change
  };
  uint8_t condense_cigar_flag[5][8] = {
      {BAM_CSOFT_CLIP, BAM_CMATCH},
      {BAM_CSOFT_CLIP, BAM_CDEL, BAM_CPAD, BAM_CDEL,
       BAM_CMATCH},                                      // change this
      {BAM_CMATCH, BAM_CPAD, BAM_CDEL, BAM_CSOFT_CLIP},  // change this
      {BAM_CMATCH, BAM_CSOFT_CLIP},
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CDEL, BAM_CMATCH, BAM_CPAD, BAM_CMATCH,
       BAM_CPAD, BAM_CMATCH}};
  uint32_t condense_cigar_len[5][8] = {{12, 139},
                                       {33, 1, 1, 1, 114},
                                       {121, 1, 4, 23},
                                       {103, 49},
                                       {23, 18, 3, 57, 1, 10, 1, 39}};
  unsigned int overlapping_primer_sizes[] = {0, 2, 2, 0, 0, 0, 0, 2, 2, 1};
  int ctr = 0;
  std::vector<primer> overlapping_primers;
  primer cand_primer;
  bool isize_flag = false;
  while (sam_itr_next(in, iter, aln) >= 0) {
    if ((aln->core.flag & BAM_FUNMAP) != 0) {
      continue;
    }
    isize_flag =
        (abs(aln->core.isize) - max_primer_len) > abs(aln->core.l_qseq);
    std::cout << bam_get_qname(aln) << std::endl;
    get_overlapping_primers(aln, primers, overlapping_primers);
    if (overlapping_primers.size() != overlapping_primer_sizes[ctr]) {
      success = -1;
      std::cout << "Overlapping primer sizes for " << bam_get_qname(aln)
                << ". Expected: " << overlapping_primer_sizes[ctr]
                << ". Got: " << overlapping_primers.size() << std::endl;
    }
    if (overlapping_primers.size() > 0) {
      print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
      if (bam_is_rev(aln)) {
        cand_primer = get_min_start(overlapping_primers);
        t = primer_trim(aln, isize_flag, cand_primer.get_start() - 1, false);
      } else {
        cand_primer = get_max_end(overlapping_primers);
        t = primer_trim(aln, isize_flag, cand_primer.get_end() + 1, false);
        aln->core.pos += t.start_pos;
      }
      if (cand_primer.get_indice() != primer_indices[primer_ctr]) {
        success = -1;
        std::cout << "Primer indice wrong. Expected: "
                  << primer_indices[primer_ctr]
                  << ". Got: " << cand_primer.get_indice() << std::endl;
      }
      // Check start pos
      if (aln->core.pos != (int)read_start_pos[primer_ctr]) {
        success = -1;
        std::cout << "Incorrect start position" << std::endl;
      }
      // Replace cigar
      replace_cigar(aln, t.nlength, t.cigar);
      cigar = bam_get_cigar(aln);
      for (uint i = 0; i < t.nlength; ++i) {
        if (((cigar[i]) & BAM_CIGAR_MASK) != cigar_flag[primer_ctr][i]) {
          std::cout << bam_cigar_op(cigar[i]) << " "
                    << bam_cigar_oplen(cigar[i]) << " i " << i << "\n";
          success = -1;
          std::cout << "Cigar flag didn't match for "
                    << cand_primer.get_indice() << " ! Expected "
                    << (uint)cigar_flag[primer_ctr][i] << " "
                    << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
        }
        if ((((cigar[i]) >> BAM_CIGAR_SHIFT)) != cigar_len[primer_ctr][i]) {
          success = -1;
          std::cout << "Cigar length didn't match for " << bam_get_qname(aln)
                    << " ! Expected " << (uint)cigar_len[primer_ctr][i] << " "
                    << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
        }
      }

      // success += flag;
      return (num_tests == success) ? 0 : -1;

      // clear memory and close files
      bam_destroy1(aln);
      bam_itr_destroy(iter);
      sam_hdr_destroy(header);
      hts_idx_destroy(idx);
      hts_close(in);
    }
