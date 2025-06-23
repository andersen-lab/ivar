#include <iostream>
#include <vector>

#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/trim_primer_quality.h"

int test_amplicon_search(std::string bam_file, std::string bed_file,
                         std::string pair_info_file, int32_t primer_offset,
                         bool expected[]) {
  std::vector<primer> primers;
  IntervalTree amplicons;
  primers = populate_from_file(bed_file, primer_offset);
  std::cout << "Total Number of primers: " << primers.size() << std::endl;
  amplicons = populate_amplicons(pair_info_file, primers);
  samFile *in = hts_open(bam_file.c_str(), "r");
  sam_hdr_t *hdr = sam_hdr_read(in);
  bam1_t *aln = bam_init1();
  int cnt = 0;
  int res = 0;
  while (sam_read1(in, hdr, aln) >= 0) {
    res = amplicon_filter(amplicons, aln);
    if (res != expected[cnt]) {
      std::cout << "Amplicon filter output did not match " << bam_get_qname(aln)
                << " Expected " << expected[cnt] << " "
                << "Got " << res << std::endl;
      std::cout << "Read intervals: " << aln->core.pos << ":" << bam_endpos(aln)
                << std::endl;
      return -1;
    }
    cnt++;
  }
  return 0;
}

int main() {
  int success;
  std::string bam = "../data/test_amplicon.sorted.bam";
  std::string pair_indice = "../data/pair_info_2.tsv";
  int32_t primer_offset = 0;
  std::string bed = "../data/test_isize.bed";
  std::string prefix = "/tmp/data/trim_amplicon";
  std::string bam_out = "/tmp/trim_amplicon.bam";

  IntervalTree amplicons;
  bool expected[8] = {true, false, true, false, true, false, true, false};
  success =
      test_amplicon_search(bam, bed, pair_indice, primer_offset, expected);
  amplicons.inOrder();
  return success;
}
