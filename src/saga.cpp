#include "saga.h"
#include "ref_seq.h"
#include "parse_gff.h"
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include "variant_caller.h"

//first main function call
int preprocess_reads(std::string bam, std::string bed, std::string bam_out, std::string cmd, std::string pair_info, int32_t primer_offset, uint32_t min_depth, uint8_t min_qual, std::string ref_file, std::string gff_path){
  if(ref_file.empty()){
    std::cerr << "Please provide reference." << std::endl;
    exit(1);
  }
  //load the reference
  ref_antd refantd(ref_file, "");

  int retval = 0;
  std::vector<primer> primers;
  if (!bed.empty()) {
    primers = populate_from_file(bed, primer_offset);
  }
  IntervalTree amplicons;
  if (!pair_info.empty() && !primers.empty()) {
    amplicons = populate_amplicons(pair_info, primers);
    amplicons.inOrder();
    amplicons.get_max_pos();
    std::cerr << "Maximum position " << amplicons.max_pos << std::endl;
  } else{
    std::cerr << "Amplicon specific measurements will not be calculated." << std::endl;
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
  // Get the header
  sam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL) {
    std::cerr << "Unable to read header from input file." << std::endl;
    return -1;
  }
  // Initiate the alignment record
  bam1_t *aln = bam_init1();
  cigar_ t;
  init_cigar(&t);
  std::vector<primer>::iterator cit;
  std::vector<bam1_t *> alns;

  //iterate through reads
  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  aln = bam_init1();

  uint32_t last_position=0;
  int tid=0;
  while (sam_read1(in, header, aln) >= 0) {
    uint32_t end_pos = aln->core.pos + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));
    if (end_pos > last_position) {
      last_position = end_pos;
      //get the region
      if(aln->core.tid > -1){
        tid = aln->core.tid;
      }
    }
  }

  const std::string ref_name = (std::string)header->target_name[tid];
  bam_destroy1(aln);
  bam_hdr_destroy(header);

  //Initiate variants
  variant_caller vc(min_qual, ref_file, gff_path);
  if(!vc.initialize_region(ref_name)) {
    std::cerr << "Unable to initialize region " << ref_name << std::endl;
    return 0;
  }
  vc.set_amplicons(amplicons);

  sam_close(in);
  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  aln = bam_init1();
  //hold the reads until it's mate can be found
  std::unordered_map<std::string, bam1_t*> read_map;
  // Iiterate through reads
  while (sam_read1(in, header, aln) >= 0) {
    //get the name of the read
    std::string read_name = bam_get_qname(aln);
    if (!(aln->core.flag & BAM_FPAIRED) || !(aln->core.flag & BAM_FPROPER_PAIR)){
      // Unpaired reads
      std::vector<site_state> read_site_states;
      vc.parse_single_read(aln, ref_name, read_site_states);
      vc.assign_amplicon_to_read(aln->core.pos, bam_endpos(aln), read_site_states);
      vc.add_variants(read_site_states);
      continue;
    }
    auto it = read_map.find(read_name);
    //assumption is that read pairs share a name
    //execute if we've already seen the mate
    if (it != read_map.end()) {
      bam1_t* mate = it->second;
      std::vector<site_state> read_site_states;
      vc.parse_paired_reads(aln, mate, ref_name, read_site_states);
      // Take care of cases where reads are out of order
      uint32_t insert_start = aln->core.pos;
      uint32_t insert_end = bam_endpos(mate);
      if(insert_start > insert_end){
        insert_start = mate->core.pos;
        insert_end = bam_endpos(aln);
      }
      vc.assign_amplicon_to_read(insert_start, insert_end, read_site_states);
      vc.add_variants(read_site_states);
      //clean the mate
      bam_destroy1(mate);
      read_map.erase(it);
    } else {
      // Store the current read in the map
      read_map[read_name] = bam_dup1(aln);  // Duplicate the read to avoid overwriting
    }
  }
  vc.write_to_file(bam_out, ref_name);
  vc.write_codon_to_file(bam_out, ref_name);
  vc.write_aa_to_file(bam_out, ref_name);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return(retval);
}
