#include <iostream>

#include "../src/trim_primer_quality.h"

int main() {
  int success = 0;
  std::string bam = "../data/primer_only/primer_edge_cases.bam";
  int32_t primer_offset = 0;
  std::vector<primer> primers =
      populate_from_file("../data/test.bed", primer_offset);
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
  cigar_ t;
  uint32_t *cigar;
  int primer_ctr = 0;
  int ctr = 0;
  int len = 0;
  std::vector<primer> overlapping_primers;
  primer cand_primer;
  bool isize_flag = false;
  while (sam_itr_next(in, iter, aln) >= 0) {
    if ((aln->core.flag & BAM_FUNMAP) != 0) {
      continue;
    }
    isize_flag =
        (abs(aln->core.isize) - max_primer_len) > abs(aln->core.l_qseq);
    len = bam_cigar2qlen(aln->core.n_cigar, bam_get_cigar(aln));
    std::cout << bam_get_qname(aln) << std::endl;
    print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    get_overlapping_primers(aln, primers, overlapping_primers, false);
    if (overlapping_primers.size() > 0) {
      // Forward trim
      cand_primer = get_max_end(overlapping_primers);
      t = primer_trim(aln, isize_flag, cand_primer.get_end() + 1, false);
      aln->core.pos += t.start_pos;
      // Replace cigar
      replace_cigar(aln, t.nlength, t.cigar);
    }
    std::cout << "Forward trim" << std::endl;
    cigar = bam_get_cigar(aln);
    print_cigar(cigar, aln->core.n_cigar);
    // Reverse trim
    get_overlapping_primers(aln, primers, overlapping_primers, true);
    if (overlapping_primers.size() > 0) {
      cand_primer = get_min_start(overlapping_primers);
      t = primer_trim(aln, isize_flag, cand_primer.get_start() - 1, true);
      aln->core.pos += t.start_pos;
      // Replace cigar
      replace_cigar(aln, t.nlength, t.cigar);
    }
    std::cout << "Reverse trim" << std::endl;
    cigar = bam_get_cigar(aln);
    print_cigar(cigar, aln->core.n_cigar);
    // Condense cigar
    std::cout << std::endl << "Condensing cigar ... " << std::endl;
    condense_cigar(&t);
    replace_cigar(aln, t.nlength, t.cigar);
    cigar = bam_get_cigar(aln);
    print_cigar(cigar, aln->core.n_cigar);
    if (bam_cigar2qlen(aln->core.n_cigar, bam_get_cigar(aln)) != len) {
      success = -1;
      std::cout << "Cigar length and read length don't match after trimming"
                << std::endl;
      std::cout << "Expected" << len << ". Got "
                << bam_cigar2qlen(aln->core.n_cigar, bam_get_cigar(aln))
                << std::endl;
    }
    if (aln->core.n_cigar != 1 || bam_cigar_op(cigar[0]) != 4) {
      success = -1;
      std::cout << "Complete primer not trimmed " << bam_cigar_op(cigar[0])
                << std::endl;
    }
    primer_ctr++;
    ctr++;
    std::cout << " ---- " << std::endl;
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
