#include <stdint.h>
#include <string.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "interval_tree.h"
#include "primer_bed.h"

#ifndef trim_primer_quality
#define trim_primer_quality

struct cigar_ {
  uint32_t *cigar;
  bool free_cig;
  uint32_t nlength;
  int32_t start_pos;
};

inline void init_cigar(cigar_ *t) {
  t->cigar = NULL;
  t->free_cig = false;
  t->nlength = 0;
  t->start_pos = 0;
}
inline void free_cigar(cigar_ t) {
  if (t.free_cig) free(t.cigar);
}

void add_pg_line_to_header(bam_hdr_t **hdr, char *cmd);

int trim_bam_qual_primer(std::string bam, std::string bed, std::string bam_out,
                         uint8_t min_qual, uint8_t sliding_window,
                         std::string cmd, bool write_no_primer_reads,
                         bool mark_qcfail_flag, int min_length,
                         std::string pair_info, int32_t primer_offset);
void free_cigar(cigar_ t);
int32_t get_pos_on_query(uint32_t *cigar, uint32_t ncigar, int32_t pos,
                         int32_t ref_start);
int32_t get_pos_on_reference(uint32_t *cigar, uint32_t ncigar, uint32_t pos,
                             uint32_t ref_start);
void reverse_qual(uint8_t *q, int l);
void reverse_cigar(uint32_t *cigar, int l);
double mean_quality(uint8_t *a, int s, int e);
cigar_ quality_trim(bam1_t *r, uint8_t qual_threshold, uint8_t sliding_window);
void print_cigar(uint32_t *cigar, int nlength);
cigar_ primer_trim(bam1_t *r, bool &isize_flag, int32_t new_pos,
                   bool unpaired_rev);
void replace_cigar(bam1_t *b, uint32_t n, uint32_t *cigar);
void condense_cigar(cigar_ *t);
void get_overlapping_primers(bam1_t *r, std::vector<primer> &primers,
                             std::vector<primer> &overlapping_primers);
void get_overlapping_primers(bam1_t *r, std::vector<primer> primers,
                             std::vector<primer> &overlapping_primers,
                             bool unpaired_rev);
int get_bigger_primer(std::vector<primer> primers);
bool amplicon_filter(IntervalTree amplicons, bam1_t *r);
std::vector<primer> insertionSort(std::vector<primer> primers, uint32_t n);
int binary_search(std::vector<primer> &primers, uint32_t target_pos, int low,
                  int high);
#endif
