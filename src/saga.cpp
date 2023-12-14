#include "saga.h"
#include "trim_primer_quality.h"
#include <fstream>
#include <cmath>

float calculate_standard_deviation(std::vector<float> frequencies) {
  float sum = 0.0, mean = 0.0, sd = 0.0;
  uint32_t i = 0;
  for(i = 0; i < frequencies.size(); ++i) {
    sum += frequencies[i];
  }
  mean = sum / frequencies.size();

  for(i = 0; i < frequencies.size(); ++i) {
    sd += pow(frequencies[i] - mean, 2);
  }
  return sqrt(sd / frequencies.size());
}

int preprocess_reads(std::string bam, std::string bed, std::string bam_out,
                             uint8_t min_qual,
                             std::string cmd,
                             std::string pair_info, int32_t primer_offset){
  /*
   * Here we're iterating every single read, and adding to the haplotype object in the amplicons schema.
   */

  std::cerr << min_qual << std::endl; //TODO handle the quality of bases
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

  // Get the header
  sam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL) {
    std::cerr << "Unable to read header from input file." << std::endl;
    return -1;
  }
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  amplicons.get_max_pos();
  std::cerr << "max pos " << amplicons.max_pos << std::endl;
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
    //TESTLINES
    if(start_pos > 3000){
      continue;
    }
    overlapping_primers.clear();
    //TODO handle unpaired
    if(strand == '+'){
      for(uint32_t i=start_pos-10; i < start_pos+10; i++){
        //std::cerr << "i " << i << " " << start_pos << std::endl;
        if (i < 0 || i > amplicons.max_pos) continue; 
        if (primer_map_forward.find(i) != primer_map_forward.end()) {
          overlapping_primers = primer_map_forward[i];
        }
        if (overlapping_primers.size() > 0) break;
      }
    }else{
      for (uint32_t i=start_pos-10; i < start_pos+10; i++)
      {
        if (i < 0 || i > amplicons.max_pos) continue;
        if (primer_map_reverse.find(i) != primer_map_reverse.end()){
          overlapping_primers = primer_map_reverse[i];
        }
        if (overlapping_primers.size() > 0) break;
      }
    }
    bam1_t *r = aln;
    //get the md tag
    uint8_t *aux = bam_aux_get(aln, "MD");
    //get the sequence
    uint8_t *seq = bam_get_seq(aln);
    //get cigar for the read
    uint32_t *cigar = bam_get_cigar(r); 
    uint32_t nlength = r->core.n_cigar; 
    uint8_t *qualities = bam_get_qual(r);
    //assign to a primer not an amplicon, because here direction matters
    //TODO handle the case of unpaired reads
    if (overlapping_primers.size() == 0){
      continue;
    }

    for(uint32_t i=0; i < overlapping_primers.size(); i++){
      uint32_t start = overlapping_primers[i].get_start();
      uint32_t end = overlapping_primers[i].get_end();

      for(uint32_t j=0; j < primers.size(); j++){
        uint32_t pstart = primers[j].get_start();
        uint32_t pend = primers[j].get_end();
        
        if (start == pstart && end == pend){
          primers[j].add_cigarotype(cigar, aln->core.pos, nlength, seq, aux, bam_get_qname(aln), qualities);
        }
      }
    }   
  }
  std::cerr << "transforming mutations" << std::endl;
  //PRIMER METHOD calculate mutations from unique cigars per primer, outputing variant frequencies
  for(uint32_t i=0; i < primers.size(); i++){
    primers[i].transform_mutations();
  }
  //exit(1);
  std::cerr << "setting amplicon level haplotypes" << std::endl;
  //AMPLICON METHOD translate this into amplicon haplotype obj of mutations per primer (ie. variant freq per amplicon)
  for (uint32_t i=0; i < primers.size(); i++){
    amplicons.set_haplotypes(primers[i]);      
  }
  std::vector<uint32_t> flagged_positions;
  
  std::cerr << "detecting variant abberations" << std::endl;
  //detect fluctuating variants
  for(uint32_t i=0; i < amplicons.max_pos; i++){
    amplicons.test_flux.clear();
    //this bit pushes all amp position vectors back to test_flux object
    amplicons.detect_abberations(i);
    if (amplicons.test_flux.size() < 2) continue;
    std::map<std::string, std::vector<float>> allele_maps;
    for(uint32_t j=0; j < amplicons.test_flux.size(); j++){
      uint32_t total_depth = amplicons.test_flux[j].depth;
      //std::cerr << "total depth " << total_depth << std::endl;
      if(total_depth < 20){
        break;
      }
      std::vector<allele> ad  = amplicons.test_flux[j].alleles;
      //std::cerr << "" << std::endl;
      //print_allele_depths(ad);
      for(uint32_t k=0; k < ad.size(); k++){
        std::string nuc = ad[k].nuc;
        uint32_t ad_depth = ad[k].depth;
        float t = (float)ad_depth / (float)total_depth; 
        if (allele_maps.find(nuc) != allele_maps.end()){
          allele_maps[nuc].push_back(t);  
        } else {
          std::vector<float> tmp;
          tmp.push_back(t);
          allele_maps[nuc] = tmp;
        }
      }
    }
    std::map<std::string, std::vector<float>>::iterator it;
    //std::cerr << "position is " << i << std::endl;
    for (it = allele_maps.begin(); it != allele_maps.end(); it++){
      /*std::cerr << it->first << std::endl;
      for(uint32_t x=0; x < it->second.size(); x++){
        std::cerr << it->second[x] << std::endl;
      }*/
      float sd = calculate_standard_deviation(it->second);
      //std::cerr << "sdv dev " << sd << std::endl;
      //TODO this is hard coded, consider it
      if (sd >= 0.05){
        flagged_positions.push_back(i);
        break;
      }
    }
    //std::cerr << "here!" << std::endl;
    //exit(1);
  }
  /*for(uint32_t i = 0; i < flagged_positions.size();i++){
    std::cerr << flagged_positions[i] << std::endl;
  }*/
  //combine amplicon counts to get total variants 
  amplicons.combine_haplotypes();
  std::vector<position> variants = amplicons.variants;
  std::cerr << "variants size " << variants.size() << std::endl;

  //write variants to a file
  ofstream file;
  file.open(bam_out + ".txt", ios::trunc);
  file << "POS\tALLELE\tDEPTH\tFREQ\tAVG_QUAL\tFLAGGED_POS\n";
  for(uint32_t i=0; i < variants.size(); i++){
    for(uint32_t j=0; j < variants[i].alleles.size(); j++){
      float freq = (float)variants[i].alleles[j].depth / (float)variants[i].depth;
      if(freq == 0){
        continue;
      }
      file << std::to_string(variants[i].pos) << "\t";
      if(variants[i].pos == 208){
        std::cerr << variants[i].alleles[j].nuc << std::endl;
      }
      file << variants[i].alleles[j].nuc << "\t";
      file << std::to_string(variants[i].alleles[j].depth) << "\t";   
      file << std::to_string(freq) << "\t";    
      file << std::to_string((float)variants[i].alleles[j].mean_qual / (float)variants[i].alleles[j].depth) << "\t";
      std::vector<uint32_t>::iterator it; 
      it = find(flagged_positions.begin(), flagged_positions.end(), variants[i].pos);
      if (it != flagged_positions.end()){
        file << "TRUE" << "\t";
      } else {
        file << "FALSE" << "\t";
      }
      file << "\n";
    }
    
  }  
  file.close();
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  //end, data has been appropriately preprocessed and problematic positions have been flagged
  //room for extension to calcualte physical linkage in the future
  return(retval);
}
