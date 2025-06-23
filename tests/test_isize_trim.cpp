#include <iostream>

#include "../src/trim_primer_quality.h"
#include "htslib/sam.h"

/*
typedef struct {
  uint8_t cigar_flag[3];
  uint32_t cigar_flag[3];
} result_t;
*/

int test_isize_trim(uint8_t min_qual, uint8_t sliding_window,
                    bool no_write_flag, bool keep_for_reanalysis,
                    int min_length, std::string testname,
                    uint8_t cigar_flag[14][3], uint32_t cigar_len[14][3]) {
  int success = 0;
  int res;
  uint32_t *cigar;

  bam1_t *aln = bam_init1();

  std::string bam = "../data/test_isize.sorted.bam";
  std::string bed = "../data/test_isize.bed";
  std::string pair_info = "";
  int32_t primer_offset = 0;
  std::string prefix = "../data/trim_isize";
  std::string bam_out = "../data/trim_isize.bam";

  std::string region_ = "";
  std::string cmd = "@PG\tID:ivar-trim\tPN:ivar\tVN:1.0.0\tCL:ivar trim\n";

  // Test and check result
  res = trim_bam_qual_primer(
      bam, bed, prefix, min_qual, sliding_window, cmd, no_write_flag,
      keep_for_reanalysis, min_length, pair_info, primer_offset);

  if (res) {
    success = -1;
    std::cerr << testname << " failed: trim_bam_qual_primer() returned " << res
              << std::endl;
  } else {
    samFile *in = hts_open(bam_out.c_str(), "r");
    if (!in) {
      success = -1;
      std::cerr << testname << " failed: Can't open " << bam_out << std::endl;
    } else {
      sam_hdr_t *hdr = sam_hdr_read(in);
      if (!hdr) {
        success = -1;
        std::cerr << testname << " failed: Can't read header from " << bam_out
                  << std::endl;
      } else {
        int n = 0;
        while (sam_read1(in, hdr, aln) >= 0) {
          /*
          if (aln->core.isize != expected[n].isize) {
              success = -1;
              std::cerr << testname << " test failed: found isize " <<
          aln->core.isize << " at record " << n << ": expected " <<
          expected[n].isize << std::endl;
              }
          */
          cigar = bam_get_cigar(aln);
          for (uint i = 0; i < aln->core.n_cigar; ++i) {
            if (((cigar[i]) & BAM_CIGAR_MASK) != cigar_flag[n][i]) {
              success = -1;
              std::cout << "Cigar flag didn't match for " << bam_get_qname(aln)
                        << " on " << i + 1 << " op"
                        << "! Expected " << (uint)cigar_flag[n][i] << " "
                        << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
            }
            if ((((cigar[i]) >> BAM_CIGAR_SHIFT)) != cigar_len[n][i]) {
              success = -1;
              std::cout << "Cigar length didn't match for "
                        << bam_get_qname(aln) << " on " << i + 1 << " op"
                        << "! Expected " << (uint)cigar_len[n][i] << " "
                        << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT)
                        << std::endl;
            }
          }
          n++;
        }
        sam_hdr_destroy(hdr);
      }
      sam_close(in);
    }
  }
  return success;
}

int main() {
  int success = 0;
  int test_res = 0;

  uint8_t cigar_flag[14][3] = {{BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CSOFT_CLIP, BAM_CMATCH},
                               {BAM_CSOFT_CLIP, BAM_CMATCH},
                               {BAM_CSOFT_CLIP, BAM_CMATCH},
                               {BAM_CSOFT_CLIP, BAM_CMATCH},
                               {BAM_CSOFT_CLIP, BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CSOFT_CLIP, BAM_CMATCH},
                               {BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CMATCH, BAM_CSOFT_CLIP},
                               {BAM_CMATCH, BAM_CSOFT_CLIP}};

  uint32_t cigar_len[14][3] = {
      {24, 363, 24}, {24, 155}, {24, 155}, {24, 64},  {24, 64},
      {24, 363, 24}, {24, 155}, {156, 24}, {156, 24}, {156, 24},
      {64, 24},      {64, 24},  {137, 10}, {137, 11}};

  test_res = test_isize_trim(20, 4, false, false, 30, "default paramaters",
                             cigar_flag, cigar_len);
  if (test_res) success = -1;

  return success;
}
/*
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7
*/
