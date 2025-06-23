#include <stdint.h>

#include <iostream>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "primer_bed.h"
#include "trim_primer_quality.h"

#ifndef removereads_from_amplicon
#define removereads_from_amplicon

int rmv_reads_from_amplicon(std::string bam, std::string region_,
                            std::string bam_out, std::vector<std::string> amp,
                            std::string bed, std::string cmd);

#endif
