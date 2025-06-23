#include <iostream>

#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/trim_primer_quality.h"
#include "htslib/sam.h"

typedef struct {
  uint16_t flag;
  uint32_t n_cigar;
  hts_pos_t isize;
  uint8_t fail;
} result_t;

int test_trim(uint8_t min_qual, uint8_t sliding_window, bool no_write_flag,
              bool keep_for_reanalysis, int min_length, std::string testname,
              int nexpected, result_t *expected) {
  int success = 0;
  int res;

  bam1_t *aln = bam_init1();

  std::string bam = "../data/test.unmapped.sorted.bam";
  std::string bed = "../data/test.bed";
  std::string pair_info = "";
  int32_t primer_offset = 0;
  std::string prefix = "/tmp/trim";
  std::string bam_out = "/tmp/trim.bam";

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
          std::cerr << (aln->core.flag & BAM_FQCFAIL ? "FAIL" : "Pass") << "\t";
          std::cerr << aln->core.n_cigar << "\t" << aln->core.l_qseq << "\t"
                    << aln->core.mtid << "\t" << aln->core.mpos << "\t"
                    << aln->core.isize << std::endl;
          if (aln->core.isize != expected[n].isize) {
            success = -1;
            std::cerr << testname << " test failed: found isize "
                      << aln->core.isize << " at record " << n << ": expected "
                      << expected[n].isize << std::endl;
          }
          if (!(aln->core.flag & BAM_FQCFAIL) != !expected[n].fail) {
            success = -1;
            std::cerr << testname << " test failed: found "
                      << (aln->core.flag & BAM_FQCFAIL ? "FAIL" : "Pass")
                      << " at record " << n << ": expected "
                      << (expected[n].fail ? "FAIL" : "Pass") << std::endl;
          }
          n++;
        }
        if (n != nexpected) {
          success = -1;
          std::cerr << testname << " test failed: found " << n
                    << " records: expected " << nexpected << std::endl;
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

  // Default paramaters
  {
    result_t expected[5] = {{163, 2, 364, 0},
                            {163, 2, 382, 0},
                            {83, 2, -356, 0},
                            {83, 2, -382, 0},
                            {67, 2, -398, 0}};
    if (test_trim(20, 4, false, false, 30, "default paramaters", 5, expected))
      success = -1;
  }

  // -e (write missing primers)
  {
    result_t expected[9] = {
        {131, 7, 403, 0}, {163, 2, 364, 0}, {163, 2, 382, 0},
        {163, 2, 242, 0}, {83, 2, -242, 0}, {147, 2, -150, 0},
        {83, 2, -356, 0}, {83, 2, -382, 0}, {67, 2, -398, 0}};
    if (test_trim(20, 4, true, false, 30, "write missing primers", 9, expected))
      success = -1;
  }

  // -k (keep for reanalysis)
  {
    result_t expected[10] = {{643, 7, 403, 1},  {163, 2, 364, 0},
                             {163, 2, 382, 0},  {675, 2, 242, 1},
                             {595, 2, -242, 1}, {611, 2, -150, 1},
                             {659, 2, -150, 1}, {83, 2, -356, 0},
                             {83, 2, -382, 0},  {67, 2, -398, 0}};
    if (test_trim(20, 4, false, true, 30, "keep for reanalysis", 10, expected))
      success = -1;
  }

  // -e -k (all flags)
  {
    result_t expected[10] = {{643, 7, 403, 0},  {163, 2, 364, 0},
                             {163, 2, 382, 0},  {675, 2, 242, 0},
                             {595, 2, -242, 0}, {611, 2, -150, 1},
                             {659, 2, -150, 0}, {83, 2, -356, 0},
                             {83, 2, -382, 0},  {67, 2, -398, 0}};
    if (test_trim(20, 4, true, true, 30, "both flags", 10, expected))
      success = -1;
  }

  return success;
}
