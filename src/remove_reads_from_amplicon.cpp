#include "remove_reads_from_amplicon.h"

int rmv_reads_from_amplicon(std::string bam, std::string region_,
                            std::string bam_out, std::vector<std::string> amp,
                            std::string bed, std::string cmd) {
  std::vector<primer> primers = populate_from_file(bed);
  if (primers.size() == 0) {
    return 0;
  }
  bam_out += ".bam";
  std::cout << "Writing to " << bam_out << std::endl;
  if (bam.empty()) {
    std::cout << "Bam in empty" << std::endl;
    return 0;
  }
  // open BAM for reading
  samFile *in = hts_open(bam.c_str(), "r");
  BGZF *out = bgzf_open(bam_out.c_str(), "w");
  if (in == NULL) {
    std::cout << ("Unable to open BAM/SAM file.") << std::endl;
    return -1;
  }
  // Load the index
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  if (idx == NULL) {
    if (sam_index_build2(bam.c_str(), 0, 0) < 0) {
      std::cout << ("Unable to open BAM/SAM index.") << std::endl;
      return -1;
    } else {
      idx = sam_index_load(in, bam.c_str());
    }
  }
  // Get the header
  bam_hdr_t *header = sam_hdr_read(in);
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  if (bam_hdr_write(out, header) < 0) {
    std::cout << "Unable to write BAM header to path." << std::endl;
    sam_close(in);
    return -1;
  }
  if (header == NULL) {
    sam_close(in);
    std::cout << "Unable to open BAM/SAM header." << std::endl;
  }
  if (region_.empty()) {
    std::cout << "Number of references: " << header->n_targets << std::endl;
    for (int i = 0; i < header->n_targets; ++i) {
      std::cout << "Reference Name: " << header->target_name[i] << std::endl;
      std::cout << "Reference Length: " << header->target_len[i] << std::endl;
      if (i == 0) {
        region_.assign(header->target_name[i]);
      }
    }
    std::cout << "Using Region: " << region_ << std::endl;
  }
  std::string temp(header->text);
  std::string sortFlag("SO:coordinate");
  if (temp.find(sortFlag)) {
    std::cout << "Sorted By Coordinate" << std::endl;  // Sort by coordinate
  } else {
    std::cout << "Not sorted" << std::endl;
  }
  // Initialize iterator
  hts_itr_t *iter = NULL;
  // Move the iterator to the region we are interested in
  iter = sam_itr_querys(idx, header, region_.c_str());
  if (header == NULL || iter == NULL) {
    sam_close(in);
    std::cout << "Unable to iterate to region within BAM/SAM." << std::endl;
    return -1;
  }
  // Initiate the alignment record
  bam1_t *aln = bam_init1();
  int ctr = 0, rmv_ctr = 0;
  bool w;
  while (sam_itr_next(in, iter, aln) >= 0) {
    uint8_t *a = bam_aux_get(aln, "XA");
    w = true;
    if (a != 0) {
      for (std::vector<std::string>::iterator it = amp.begin(); it != amp.end();
           ++it) {
        if (bam_aux2i(a) == get_primer_indice(primers, *it)) {
          w = false;
        }
      }
    }
    if (w) {
      if (bam_write1(out, aln) < 0) {
        std::cout << "Not able to write to BAM" << std::endl;
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        bam_destroy1(aln);
        bam_hdr_destroy(header);
        sam_close(in);
        bgzf_close(out);
        return -1;
      };
    } else {
      rmv_ctr++;
    }
    ctr++;
    if (ctr % 1000000 == 0) {
      std::cout << "Processed " << ctr << " reads" << std::endl;
    }
  }
  std::cout << "Results:" << std::endl;
  std::cout << rmv_ctr << " reads were removed." << std::endl;
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  bgzf_close(out);
  return 0;
}
