#include <iostream>

#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/trim_primer_quality.h"
#include "htslib/sam.h"

int main() {
  int success = 0;
  std::string bam = "../data/test.sim.merged.sorted.bam";
  std::vector<primer> primers = populate_from_file("../data/test_merged.bed");
  int max_primer_len = 0;
  max_primer_len = get_bigger_primer(primers);
  std::string region_;
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
  bam_hdr_t *header = sam_hdr_read(in);
  region_.assign(header->target_name[0]);
  std::string temp(header->text);
  hts_itr_t *iter = NULL;
  iter = sam_itr_querys(idx, header, region_.c_str());
  bam1_t *aln = bam_init1();
  cigar_ t = {NULL, false, 0, 0};
  uint32_t *cigar;
  int primer_ctr = 0;
  int forward_primer_indices[] = {1, 1, 5, 5, 6};
  int rev_primer_indices[] = {3, -1, -1, -1, 4};
  uint8_t cigar_flag[5][17] = {
      {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CDEL, BAM_CMATCH, BAM_CDEL,
       BAM_CMATCH, BAM_CINS, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CPAD,
       BAM_CMATCH, BAM_CDEL, BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP,
       BAM_CSOFT_CLIP},  // 3S15M5D2M1D14M1I56M1P6M3P93M2D4M2I8M
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CPAD, BAM_CMATCH,
       BAM_CPAD, BAM_CMATCH, BAM_CPAD, BAM_CMATCH},
      {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP,
       BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CSOFT_CLIP},
      {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CDEL, BAM_CPAD,
       BAM_CDEL, BAM_CMATCH, BAM_CPAD, BAM_CMATCH},
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CSOFT_CLIP}};
  uint32_t cigar_len[5][17] = {
      {3, 15, 5, 2, 1, 14, 1, 56, 1, 6, 3, 93, 2, 1, 3, 2, 8},  // change here
      {1, 19, 1, 56, 1, 6, 3, 99, 2, 65},
      {8, 13, 1, 6, 3, 99, 2, 62, 3},
      {9, 10, 3, 3, 3, 2, 97, 2, 26},  // change here
      {18, 71, 2, 89, 18}};
  uint32_t read_start_pos[5] = {89, 92, 173, 170, 201};
  uint8_t condense_cigar_flag[5][14] = {
      {BAM_CSOFT_CLIP, BAM_CDEL, BAM_CMATCH, BAM_CDEL, BAM_CMATCH, BAM_CINS,
       BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CDEL,
       BAM_CMATCH, BAM_CSOFT_CLIP},
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CPAD, BAM_CMATCH,
       BAM_CPAD, BAM_CMATCH, BAM_CPAD, BAM_CMATCH},
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CSOFT_CLIP},
      {BAM_CSOFT_CLIP, BAM_CDEL, BAM_CPAD, BAM_CDEL, BAM_CMATCH, BAM_CPAD,
       BAM_CMATCH},  // change here
      {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CPAD, BAM_CMATCH, BAM_CSOFT_CLIP}};
  uint32_t condense_cigar_len[5][14] = {
      {18, 5, 2, 1, 14, 1, 56, 1, 6, 3, 93, 2, 1, 13},  // change here
      {1, 19, 1, 56, 1, 6, 3, 99, 2, 65},
      {31, 99, 2, 62, 3},
      {22, 3, 3, 2, 97, 2, 26},  // change here
      {18, 71, 2, 89, 18}};
  unsigned int fwd_overlapping_primer_sizes[] = {2, 1, 1, 1, 1};
  unsigned int rev_overlapping_primer_sizes[] = {2, 0, 0, 0, 1};
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
    std::cout << std::endl << "Forward" << std::endl;
    get_overlapping_primers(aln, primers, overlapping_primers, false);
    // Forward primer
    if (overlapping_primers.size() > 0) {
      cand_primer = get_max_end(overlapping_primers);
      t = primer_trim(aln, isize_flag, cand_primer.get_end() + 1, false);
      aln->core.pos += t.start_pos;
      replace_cigar(aln, t.nlength, t.cigar);
      if (overlapping_primers.size() !=
          fwd_overlapping_primer_sizes[primer_ctr]) {
        success = -1;
        std::cout << "Overlapping primer sizes for " << bam_get_qname(aln)
                  << ". Expected: " << fwd_overlapping_primer_sizes[primer_ctr]
                  << ". Got: " << overlapping_primers.size() << std::endl;
      }
      if (cand_primer.get_indice() != forward_primer_indices[primer_ctr]) {
        success = -1;
        std::cout << "Primer indice wrong. Expected: "
                  << forward_primer_indices[primer_ctr]
                  << ". Got: " << cand_primer.get_indice() << std::endl;
      }
    }
    // Reverse primer
    std::cout << std::endl << "Reverse" << std::endl;
    get_overlapping_primers(aln, primers, overlapping_primers, true);
    if (overlapping_primers.size() > 0) {
      cand_primer = get_min_start(overlapping_primers);
      free_cigar(t);
      t = primer_trim(aln, isize_flag, cand_primer.get_start() - 1, true);
      replace_cigar(aln, t.nlength, t.cigar);
      if (overlapping_primers.size() !=
          rev_overlapping_primer_sizes[primer_ctr]) {
        success = -1;
        std::cout << "Overlapping primer sizes for " << bam_get_qname(aln)
                  << ". Expected: " << rev_overlapping_primer_sizes[primer_ctr]
                  << ". Got: " << overlapping_primers.size() << std::endl;
      }
      if (cand_primer.get_indice() != rev_primer_indices[primer_ctr]) {
        success = -1;
        std::cout << "Primer indice wrong. Expected: "
                  << rev_primer_indices[primer_ctr]
                  << ". Got: " << cand_primer.get_indice() << std::endl;
      }
    }
    if (aln->core.pos != (int)read_start_pos[primer_ctr]) {
      success = -1;
      std::cout << "Start pos didn't match" << std::endl;
    }
    replace_cigar(aln, t.nlength, t.cigar);
    cigar = bam_get_cigar(aln);
    for (uint i = 0; i < t.nlength; ++i) {
      if (((cigar[i]) & BAM_CIGAR_MASK) != cigar_flag[primer_ctr][i]) {
        success = -1;
        std::cout << "Cigar flag didn't match for " << bam_get_qname(aln)
                  << " ! Expected " << (uint)cigar_flag[primer_ctr][i] << " "
                  << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
      }
      if ((((cigar[i]) >> BAM_CIGAR_SHIFT)) != cigar_len[primer_ctr][i]) {
        success = -1;
        std::cout << "Cigar length didn't match for " << bam_get_qname(aln)
                  << " ! Expected " << (uint)cigar_len[primer_ctr][i] << " "
                  << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
      }
    }
    // Condense cigar
    std::cout << std::endl << "Condensing cigar ... " << std::endl;
    condense_cigar(&t);
    replace_cigar(aln, t.nlength, t.cigar);
    cigar = bam_get_cigar(aln);
    for (uint i = 0; i < t.nlength; ++i) {
      if (((cigar[i]) & BAM_CIGAR_MASK) != condense_cigar_flag[primer_ctr][i]) {
        success = -1;
        std::cout << "Cigar flag didn't match! Expected "
                  << condense_cigar_flag[primer_ctr][i] << " "
                  << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
      }
      if ((((cigar[i]) >> BAM_CIGAR_SHIFT)) !=
          condense_cigar_len[primer_ctr][i]) {
        success = -1;
        std::cout << "Cigar length didn't match after condense! Expected "
                  << condense_cigar_len[primer_ctr][i] << " "
                  << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
      }
    }
    primer_ctr++;
    free_cigar(t);
  }
  // Check if primers found at all
  success = (primer_ctr > 0) ? success : -1;

  bam_destroy1(aln);
  bam_itr_destroy(iter);
  sam_hdr_destroy(header);
  hts_idx_destroy(idx);
  hts_close(in);

  return success;
}
