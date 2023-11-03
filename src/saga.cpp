#include "saga.h"
#include "trim_primer_quality.h"



int preprocess_reads(std::string bam, std::string bed, std::string bam_out,
                             uint8_t min_qual,
                             std::string cmd,
                             std::string pair_info, int32_t primer_offset){
  /*
   * Here we're iterating every single read, and adding to the haplotype object in the amplicons schema.
   */

  std::cerr << min_qual << std::endl;
  int retval = 0;
  std::vector<primer> primers;
  int max_primer_len = 0;
  if (!bed.empty()) {
    primers = populate_from_file(bed, primer_offset);
    if (primers.size() == 0) {
      std::cerr << "Exiting." << std::endl;
      return -1;
    }
  }

  // calculate the primers that should cover each position
  std::vector<std::map<uint32_t, std::vector<primer>>> hash = find_primer_per_position(primers);
  std::map<uint32_t, std::vector<primer>> primer_map_forward = hash[0];
  std::map<uint32_t, std::vector<primer>> primer_map_reverse = hash[1];

  max_primer_len = get_bigger_primer(primers);
  // get coordinates of each amplicon
  IntervalTree amplicons;
  if (!pair_info.empty()) {
    amplicons = populate_amplicons(pair_info, primers);
    /*if (amplicons.size() == 0){
      std::cerr << "Exiting." << std::endl;
      return -1;
    }*/
    std::cerr << "Amplicons detected: " << std::endl;
    amplicons.inOrder();
  }

  // Read in input file
  samFile *in;
  
  if(bam.empty()) {
    std::cerr << "Reading from stdin" << std::endl;
    in = sam_open("-", "r");
  } else {
    in = sam_open(bam.c_str(), "r");
    std::cerr << "Reading from " << bam << std::endl;
  }
  if (in == NULL) {
    std::cerr << ("Unable to open input file.") << std::endl;
    return -1;
  }

  // Setup output file
  samFile *out;
  if(bam_out.empty()) {
    out = sam_open("-", "w");
  } else {
    bam_out += ".bam";
    out = sam_open(bam_out.c_str(), "wb");
  }

  // Get the header
  sam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL) {
    std::cerr << "Unable to read header from input file." << std::endl;
    return -1;
  }
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  if (sam_hdr_write(out, header) < 0) {
    std::cerr << "Unable to write BAM header to path." << std::endl;
    sam_close(in);
    return -1;
  }

  // Initiate the alignment record
  bam1_t *aln = bam_init1();
  //int ctr = 0;
  cigar_ t;
  init_cigar(&t);
  //bool unmapped_flag = false;
  //bool amplicon_flag = false;
  //primer cand_primer;
  //bool isize_flag = true;
  //uint32_t failed_frag_size = 0;
  //uint32_t unmapped_counter = 0;
  //uint32_t amplicon_flag_ctr = 0;
  std::vector<primer>::iterator cit;
  std::vector<primer> overlapping_primers;
  std::vector<bam1_t *> alns;

  //int cig;
  char strand = '+';
  uint32_t start_pos = -1;

  // Iterate through reads
  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  aln = bam_init1();
 
  // Iterate through reads
  while (sam_read1(in, header, aln) >= 0) {
    strand = '+';
    if (bam_is_rev(aln)) {
      start_pos = bam_endpos(aln) - 1;
      strand = '-';
    } else {
      start_pos = aln->core.pos;
    }
    overlapping_primers.clear();
    //for this case, we've already trimmed so the starting pos will be shifted
    //TODO look instead for primers matching a small range!
    if(strand == '+'){
      if (primer_map_forward.find(start_pos) != primer_map_forward.end()) {
        overlapping_primers = primer_map_forward[start_pos];
      }
    }else{
      if (primer_map_reverse.find(start_pos) != primer_map_reverse.end()){
        overlapping_primers = primer_map_reverse[start_pos];
      }
    }
    //get cigar for the read
    bam1_t *r = aln;
    uint32_t *cigar = bam_get_cigar(r); 
    uint32_t nlength = r->core.n_cigar; 
    
    //assign to a primer not an amplicon, because here direction matters
    if (overlapping_primers.size() == 0){
      //std::cerr << start_pos << " " << strand << std::endl;
      //std::cerr << "this" << std::endl;
      continue;
    }

    for(uint32_t i=0; i < overlapping_primers.size(); i++){
      uint32_t start = overlapping_primers[i].get_start();
      uint32_t end = overlapping_primers[i].get_end();

      for(uint32_t j=0; j < primers.size(); j++){
        uint32_t pstart = primers[j].get_start();
        uint32_t pend = primers[j].get_end();

        if (start == pstart && end == pend){
          primers[j].add_cigarotype(cigar, start_pos, nlength);
        }
      }
    }   
  }

  for(uint32_t i=0; i < primers.size(); i++){
    std::cout << "primer " << i << std::endl;
    primer tmp = primers[i];
    std::vector<std::vector<uint32_t>> otmp = tmp.get_cigarotypes();
    std::vector<uint32_t> nlengths = tmp.get_nlengths();
    for(uint32_t j=0; j < otmp.size(); j++){
      std::cerr << "cigarotype " << j << std::endl;
      std::cerr << nlengths[j] << std::endl;

      for(uint32_t k = 0; k < nlengths[j]; k++){
        std::cerr << "k " << k << std::endl;
        std:: cerr << bam_cigar_op(otmp[j][k]) << " " << bam_cigar_oplen(otmp[j][k]) << std::endl;
      }
    }
    if( i > 10) break;
  }
  
  //PRIMER METHOD calculate mutations from unique cigars per primer (ie. count mutations per primer efficiently)
  //AMPLICON METHOD translate this into amplicon haplotype obj of mutations per primer (ie. variant freq per amplicon)
  //detect fluctuating variants - iterate every position and look for fluctuation between every amplicon objects, flag these
  //combine amplicon counts to get total variants
  //end, data has been appropriately preprocessed and problematic positions have been flagged
  //room for extension to calcualte physical linkage in the future
  return(retval);
}
