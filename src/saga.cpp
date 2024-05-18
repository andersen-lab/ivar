#include "saga.h"
#include "ref_seq.h"
#include "trim_primer_quality.h"
#include <fstream>
#include <cmath>
#include <chrono>
#include <tuple>
using namespace std::chrono;

void generate_range(uint32_t start, uint32_t end, std::vector<uint32_t> &result) {
  if (start > end) {
    std::swap(start, end);
  }
  for (uint32_t i = start; i <= end; ++i) {
    result.push_back(i);
  }
}


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

//first main function call
int preprocess_reads(std::string bam, std::string bed, std::string bam_out,
                             std::string cmd,
                             std::string pair_info, int32_t primer_offset){
 
  int retval = 0;
  std::vector<primer> primers;
  if (!bed.empty()) {
    primers = populate_from_file(bed, primer_offset);
    if (primers.size() == 0) {
      std::cerr << "Exiting." << std::endl;
      return -1;
    }
  }
  //load the reference sequence
  std::string gff_path = "";
  std::string ref_path = "/home/chrissy/sequence.fasta";
  ref_antd refantd(ref_path, gff_path);  

  // calculate the primers that should cover each position
  std::vector<std::map<uint32_t, std::vector<primer>>> hash = find_primer_per_position(primers);
  std::map<uint32_t, std::vector<primer>> primer_map_forward = hash[0];
  std::map<uint32_t, std::vector<primer>> primer_map_reverse = hash[1];

  //int max_primer_len = get_bigger_primer(primers);
  // get coordinates of each amplicon
  IntervalTree amplicons;
  if (!pair_info.empty()) {
    amplicons = populate_amplicons(pair_info, primers);
    std::cerr << "Amplicons detected: " << std::endl;
    amplicons.inOrder();
    amplicons.get_max_pos();
    amplicons.populate_variants();
    std::cerr << "Maximum position " << amplicons.max_pos << std::endl;
  } else{
    std::cerr << "Exiting." << std::endl;
    return -1;
  }
  /*if(amplicons.size() == 0){
    std::cerr << "Exiting." << std::endl;
    return(-1);
  }*/

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
  // Initiate the alignment record
  bam1_t *aln = bam_init1();
  //int ctr = 0;
  cigar_ t;
  init_cigar(&t);
  //bool unmapped_flag = false;
  //bool amplicon_flag = false;
  //bool isize_flag = true;
  //uint32_t failed_frag_size = 0;
  //uint32_t unmapped_counter = 0;
  //uint32_t amplicon_flag_ctr = 0;
  std::vector<primer>::iterator cit;
  std::vector<primer> overlapping_primers;
  std::vector<bam1_t *> alns;
   
  std::vector<primer> unpaired_primers;
  for(uint32_t i=0; i < primers.size(); i++){
    bool paired = amplicons.unpaired_primers(primers[i]);
    if(!paired){
      unpaired_primers.push_back(primers[i]);
    }
  }

  char strand = '+';
  uint32_t start_pos = -1;
  uint32_t outside_amp = 0;
  uint32_t lower_search=0;
  uint32_t counter = 0;
  // Iterate through reads
  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  aln = bam_init1();
  // Iterate through reads
  std::string test = "";
  while (sam_read1(in, header, aln) >= 0) {
    strand = '+';
    if (bam_is_rev(aln)) {
      start_pos = bam_endpos(aln) - 1;
      strand = '-';
    } else {
      start_pos = aln->core.pos;
    }
    counter += 1;
    overlapping_primers.clear();
    if(strand == '+'){
      if(start_pos >= 10){
        lower_search = start_pos-10;
      } else{
        lower_search = 0;
      }
      for(uint32_t i=lower_search; i < start_pos+10; i++){
        if (i > amplicons.max_pos) continue; 
        if (primer_map_forward.find(i) != primer_map_forward.end()) {
          overlapping_primers = primer_map_forward[i];
        }
        if (overlapping_primers.size() > 0) break;
      }
    }else{
      if(start_pos >= 10){
        lower_search = start_pos-10;
      } else{
        lower_search = 0;
      }
      for (uint32_t i=lower_search; i < start_pos+10; i++){
        if (i > amplicons.max_pos) continue;
        if (primer_map_reverse.find(i) != primer_map_reverse.end()){
          overlapping_primers = primer_map_reverse[i];
        }
        if (overlapping_primers.size() > 0) break;
      }
    }
    uint32_t ref_id = aln->core.tid;
    std::string region = header->target_name[ref_id];
    bam1_t *r = aln;
    //get the md tag
    uint8_t *aux = bam_aux_get(aln, "MD");
    uint32_t k=0;
    do{
      k++;
    } while(aux[k] != '\0');
    if(k <= 1){
      return(0);
    }
    //get the sequence
    uint8_t *seq = bam_get_seq(aln);
    //get cigar for the read
    uint32_t *cigar = bam_get_cigar(r); 
    uint32_t nlength = r->core.n_cigar; 
    uint8_t *qualities = bam_get_qual(r);
        
    //this handles the case of reads outside of an amplicon
    if (overlapping_primers.size() == 0){
      amplicons.add_read_variants(cigar, aln->core.pos, nlength, seq, aux, qualities, bam_get_qname(aln));
      outside_amp += 1;
      continue;
    }
    for(uint32_t i=0; i < overlapping_primers.size(); i++){
      uint32_t start = overlapping_primers[i].get_start();
      uint32_t end = overlapping_primers[i].get_end();
      bool cont = true;
      for(uint32_t k =0; k < unpaired_primers.size(); k++){
        if (unpaired_primers[k].get_start() == start && unpaired_primers[k].get_end() == end){
          amplicons.add_read_variants(cigar, aln->core.pos, nlength, seq, aux, qualities, bam_get_qname(aln));         
          outside_amp += 1;
          cont = false;
          continue;
        }
      }
      if(!cont){
        continue;
      }
      for(uint32_t j=0; j < primers.size(); j++){
        uint32_t pstart = primers[j].get_start();
        uint32_t pend = primers[j].get_end(); 
        if (start == pstart && end == pend){
          primers[j].add_cigarotype(cigar, aln->core.pos, nlength, seq, aux, bam_get_qname(aln), qualities);
        }
      }
    }
  }
  std::cerr << "counter " << counter << std::endl;
  std::cerr << "number of reads outside an amplicon: " << outside_amp << std::endl;
  //this is super time costly

  for(uint32_t i=0; i < primers.size(); i++){
    primers[i].populate_positions();
    primers[i].transform_mutations(ref_path);
  }
  std::cerr << "setting amplicon level haplotypes" << std::endl;
  for (uint32_t i=0; i < primers.size(); i++){
    amplicons.set_haplotypes(primers[i]);      
  }
  //exit(1);  
  //combine amplicon counts to get total variants 
  amplicons.combine_haplotypes();
  //detect primer binding issues
  std::vector<position> variants = amplicons.variants;
  //add in primer info
  for(uint32_t i=0; i < variants.size(); i++){
    bool mutation = false;
    //establish a mutation at this position
    for(auto al : variants[i].alleles){
      if(!al.is_ref && (((float)al.depth/(float)variants[i].depth) > 0.05)){
        mutation = true;
      }
    }
    if(!mutation) continue;
    //establish the mutation is within a primer region
    for(uint32_t j=0; j < primers.size(); j++){
      if(variants[i].pos >= primers[j].get_start() && variants[i].pos <= primers[j].get_end()+1){
        if(primers[j].get_strand() == '+'){
          //find the amplicon associated with this primer
          amplicons.detect_primer_issues(primers[j].get_start());
        } else {
          amplicons.detect_primer_issues(primers[j].get_end() + 1);
        }
      }
    }
  }
  //then when we experience flux on primer bound issue amplicon, save adjusted depths
  std::vector<uint32_t> flagged_positions; 
  std::vector<std::tuple<uint32_t, uint32_t, std::string>> adjusted_depths; //holds adjusted depth info
  std::cerr << "detecting variant abberations" << std::endl;
  uint32_t test_pos = 0;
  //detect fluctuating variants
  for(uint32_t i=0; i < amplicons.max_pos; i++){
    amplicons.test_flux.clear();
    amplicons.test_test.clear();
    
    //this bit pushes all amp position vectors back to test_flux object
    amplicons.detect_abberations(i);
    if(i == test_pos){
      for(auto xx : amplicons.test_test){
        std::cerr << "test flux " << xx << std::endl;
      }
    }
    if (amplicons.test_flux.size() < 2) continue;
    
    std::map<std::string, std::vector<float>> allele_maps;
    for(uint32_t j=0; j < amplicons.test_flux.size(); j++){
      if(i == test_pos){
        std::cerr << "\n" << amplicons.test_test[j] << std::endl;
      }
      uint32_t total_depth = amplicons.test_flux[j].depth;
      if(total_depth < 20){
         continue;
      }
      std::vector<allele> ad  = amplicons.test_flux[j].alleles;
      for(uint32_t k=0; k < ad.size(); k++){
        std::string nuc = ad[k].nuc;
        uint32_t ad_depth = ad[k].depth;
        if(i == test_pos){
          std::cerr << nuc << " " << ad_depth << " " << total_depth << std::endl;
        }
        float t = (float)ad_depth / (float)total_depth; 
        if (allele_maps.find(nuc) != allele_maps.end()){
          if(i == test_pos){
            std::cerr << "allele map " << nuc << " " << t << std::endl;
           }
          allele_maps[nuc].push_back(t);  
        } else {
          std::vector<float> tmp;
          tmp.push_back(t);
          if(i == test_pos){
            std::cerr << "allele map " << nuc << " " << t << std::endl;
          }
          allele_maps[nuc] = tmp;
        }
      }
    }
    std::map<std::string, std::vector<float>>::iterator it;
    for (it = allele_maps.begin(); it != allele_maps.end(); it++){
      float sd = calculate_standard_deviation(it->second);
      if(i == test_pos){
        std::cerr << i << " " << sd << " " << it->first << std::endl;
      }
      //TODO this is hard coded, consider it
      if (sd >= 0.03){
        for(auto overlap : amplicons.overlaps){
            //identify amplicons with primer binding for this site
            if(i < overlap[1] && i > overlap[0]){
              //index this amplicon within the allele map
              auto it = std::find(amplicons.test_test.begin(), amplicons.test_test.end(), overlap[0]);
              uint32_t index;
              if (it != amplicons.test_test.end()) {
                index = (uint32_t)std::distance(amplicons.test_test.begin(), it);
                position removal = amplicons.test_flux[index];
                for(auto ad : removal.alleles){
                  adjusted_depths.push_back(std::make_tuple(i, ad.depth, ad.nuc));
                  //std::cerr << "remove this depth " << ad.nuc << " " << ad.depth << std::endl;
                }
              }
            }
        }
        flagged_positions.push_back(i);
        break;
      }
    }
  }
  //exit(1);
  //TESTLINES print out the adjusted depth for sanity check
  /*
  for(uint32_t k=0; k < adjusted_depths.size(); k++){
    uint32_t pos = std::get<0>(adjusted_depths[k]);
    uint32_t depth = std::get<1>(adjusted_depths[k]);
    std::string nuc = std::get<2>(adjusted_depths[k]);
    std::cerr << pos << " " << nuc << " " << depth << std::endl;
  }*/
  //exit(1);
 
  std::vector<uint32_t> primer_binding_error;
  for(uint32_t i=0; i < amplicons.overlaps.size(); i++){
    generate_range(amplicons.overlaps[i][0], amplicons.overlaps[i][1], primer_binding_error);    
  } 

  std::cerr << "variants size " << variants.size() << std::endl;
  amplicons.overlaps.clear();
  //find amplicons that are problematic
  for(uint32_t i=0; i < flagged_positions.size(); i++){
    amplicons.detect_amplicon_overlaps(flagged_positions[i]);
  }
  std::vector<uint32_t> amplicon_level_error;
  for(uint32_t i=0; i < amplicons.overlaps.size(); i++){
    generate_range(amplicons.overlaps[i][0], amplicons.overlaps[i][1], amplicon_level_error);    
  } 

  //write variants to a file
  ofstream file;
  file.open(bam_out + ".txt", ios::trunc);
  file << "POS\tALLELE\tDEPTH\tFREQ\tAVG_QUAL\tFLAGGED_POS\tAMP_MASKED\tPRIMER_BINDING\tADJUSTED_DEPTH\tREF\n";
  for(uint32_t i=0; i < variants.size(); i++){
    for(uint32_t j=0; j < variants[i].alleles.size(); j++){
      uint32_t adjusted_depth = 0;
      for(uint32_t k=0; k < adjusted_depths.size(); k++){
          uint32_t pos = std::get<0>(adjusted_depths[k]);
          uint32_t depth = std::get<1>(adjusted_depths[k]);
          std::string nuc = std::get<2>(adjusted_depths[k]);
          if(pos == variants[i].pos && nuc == variants[i].alleles[j].nuc){
            adjusted_depth = depth;
          }
      } 
      float freq = (float)variants[i].alleles[j].depth / (float)variants[i].depth;
      if((float)variants[i].alleles[j].depth == 0){
        continue;
      }
      file << std::to_string(variants[i].pos) << "\t";
      std::string test_string = variants[i].alleles[j].nuc;
      file << variants[i].alleles[j].nuc << "\t";
      file << std::to_string(variants[i].alleles[j].depth) << "\t";   
      file << std::to_string(freq) << "\t";    
      file << std::to_string((float)variants[i].alleles[j].mean_qual / (float)variants[i].alleles[j].depth) << "\t";
      std::vector<uint32_t>::iterator it; 
      it = find(flagged_positions.begin(), flagged_positions.end(), variants[i].pos);
      if (it != flagged_positions.end()){
        file << "TRUE\t";
      } else {
        file << "FALSE\t";
      }
      it = find(amplicon_level_error.begin(), amplicon_level_error.end(), variants[i].pos);
      if(it != amplicon_level_error.end()){
        file << "TRUE\t";
      } else {
        file << "FALSE\t";
      }
      it = find(primer_binding_error.begin(), primer_binding_error.end(), variants[i].pos);
      if(it != primer_binding_error.end()){
        file << "TRUE\t";
      } else {
        file << "FALSE\t";
      }
      //insert the adjusted depth into the file for potential primer binding locations
      file << std::to_string(variants[i].alleles[j].depth-adjusted_depth) << "\t";
      if (variants[i].alleles[j].is_ref){
        file << "TRUE";
      } else {
        file << "FALSE";
      }
      file << "\n";
    } 
  }  
  file.close();
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return(retval);
}
