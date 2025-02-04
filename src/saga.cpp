#include "saga.h"
#include "trim_primer_quality.h"
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <chrono>
#include <tuple>
using namespace std::chrono;

void parse_cigar(const bam1_t* read1, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint8_t> &qualities, uint32_t total_ref_pos, uint32_t total_query_pos, uint32_t ref_start_pos, std::vector<std::string> &overlap_seq, uint32_t o_lower, uint32_t o_upper){
  //ref_start_pos describe the point after which we start recording bases
  const uint8_t* seq_field1 = bam_get_seq(read1);
  uint32_t *cigar1 = bam_get_cigar(read1);
  uint8_t* qual = bam_get_qual(read1);

  for (uint32_t i = 0; i < read1->core.n_cigar; i++){
    uint32_t op = bam_cigar_op(cigar1[i]);
    uint32_t len = bam_cigar_oplen(cigar1[i]);    
    uint32_t counter = 0;
    
    if(op == 0){ 
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        std::string tmp = std::string(1, nuc);
        if(total_ref_pos+counter >= ref_start_pos){
          positions.push_back(total_ref_pos+counter);
          bases.push_back(tmp);
          qualities.push_back(qual[j]);
          //std::cerr << total_ref_pos+counter << " " << nuc << std::endl;      
        }
        if(total_ref_pos+counter >= o_lower && total_ref_pos+counter <= o_upper){
          overlap_seq.push_back(tmp);
        }
        counter++;
      }
    }
    if(op == 2){    
       for(uint32_t j=0; j < len; j++){
        char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        if(total_ref_pos+counter >= ref_start_pos){
          positions.push_back(total_ref_pos+counter);
          bases.push_back("-");
          qualities.push_back((uint8_t)20);
          //std::cerr << total_ref_pos+counter << " - " << std::endl;      
        }
        if(total_ref_pos+counter >= o_lower && total_ref_pos+counter <= o_upper){
          overlap_seq.push_back("-");
        }
        counter++;
      }
    }
    if(op == 1){
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        std::string tmp = "+" + std::string(1, nuc);
        if(total_ref_pos+counter >= ref_start_pos){
          positions.push_back(total_ref_pos+counter);
          bases.push_back(tmp);
          qualities.push_back(qual[j]);
          //std::cerr << total_ref_pos << " " << nuc << std::endl;
        }
        if(total_ref_pos+counter >= o_lower && total_ref_pos+counter <= o_upper){
          overlap_seq.push_back(tmp);
        }
        counter++;
      }
    }

    //consumes ref
    if(bam_cigar_type(op) & 2){
      total_ref_pos += len;
    }
    if(bam_cigar_type(op) & 1){
      total_query_pos += len;
    }
  }  
}

std::vector<uint8_t> get_aux(const bam1_t* aln){
  uint8_t *aux = bam_aux_get(aln, "MD");
  std::vector<uint8_t> aux_reformat;
  uint32_t i = 0;
  do{
    aux_reformat.push_back(aux[i]);
    i++;
  } while(aux[i] != '\0');
  return aux_reformat; 
}

uint32_t find_sequence_end(const bam1_t* read){
  uint32_t start = read->core.pos;
  uint32_t *cigar = bam_get_cigar(read); 
  //find the end of the forward read
  for (uint32_t i = 0; i < read->core.n_cigar; i++){
    uint32_t op = bam_cigar_op(cigar[i]);
    uint32_t len = bam_cigar_oplen(cigar[i]);
    //consumes ref
    if(bam_cigar_type(op) & 2){
      start += len;
    }  
  }
  return start;
}

void merge_reads(const bam1_t* read1, const bam1_t* read2, IntervalTree amplicons){
  //pass the forward first then reverse 
  //underlying assumption here is that the overlap region is identical
  //also assumes that the forward read starts more "left" than  the reverse

  //get coordinates for potential overlap area
  uint32_t end_forward = find_sequence_end(read1);
  uint32_t start_reverse = read2->core.pos;
 
  uint32_t start_forward = read1->core.pos;
  uint32_t end_reverse = find_sequence_end(read2);
 
  //iterate the first cigar string
  uint32_t total_ref_pos = start_forward;
  uint32_t total_query_pos = 0;
 
  //record the positions and their bases
  std::vector<uint32_t> positions;
  std::vector<std::string> bases;
  std::vector<uint8_t> qualities;
 
  std::vector<std::string> overlap_seq1;
  std::vector<std::string> overlap_seq2; 

  //we use all of the first read 
  parse_cigar(read1, positions, bases, qualities, total_ref_pos, total_query_pos, start_forward, overlap_seq1, start_reverse, end_forward);
  parse_cigar(read2, positions, bases, qualities, start_reverse, total_query_pos, end_forward, overlap_seq2, start_reverse, end_forward);  

  //check to make sure the read overlap is the same
  /*for(uint32_t j=0; j < overlap_seq1.size(); j++){
    if(overlap_seq1[j] != overlap_seq2[j]){
      std::cerr << "HERE"<<std::endl;
      std::cerr << overlap_seq1[j] << " " << overlap_seq2[j] << std::endl;
    }
  }*/
  
  //find assigned amplicon and populate position vector
  amplicons.find_read_amplicon(start_forward, end_reverse, positions, bases, qualities);    
}

void generate_range(uint32_t start, uint32_t end, std::vector<uint32_t> &result) {
  if (start > end) {
    std::swap(start, end);
  }
  for (uint32_t i = start; i <= end; ++i) {
    result.push_back(i);
  }
}

float calculate_standard_deviation(std::vector<float> frequencies, std::vector<uint32_t> depths) {
  double weighted_sum = 0.0, total_weight = 0.0;
  for(uint32_t i = 0; i < frequencies.size(); i++) {
    weighted_sum += frequencies[i] * depths[i];
    total_weight += depths[i];
  }
  double mean = weighted_sum / total_weight;
  double weighted_variance = 0.0;
  for (uint32_t i = 0; i < frequencies.size(); i++) {
    weighted_variance += depths[i] * std::pow(frequencies[i] - mean, 2);
  }
  weighted_variance /= total_weight;
  return std::sqrt(weighted_variance);
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
  std::string gff_path = "";

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
    //amplicons.amplicon_position_pop();
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

  //hold the reads until it's mate can be found
  std::unordered_map<std::string, bam1_t*> read_map;
   

  // Iiterate through reads
  while (sam_read1(in, header, aln) >= 0) {
    strand = '+';  
    if (bam_is_rev(aln)) {
      start_pos = bam_endpos(aln) - 1;
      strand = '-';
    } else {
      start_pos = aln->core.pos;
    }
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


    /*
    if (aln->core.flag & BAM_FPAIRED){
      //std::cerr << bam_get_qname(aln) << std::endl;

    } else{
      //TODO HANDLE THIS CASE
      std::cerr << "read is unpaired" << std::endl;
      std::cerr << bam_get_qname(aln) << std::endl;
      exit(0);    
    }

    //get the name of the read    
    std::string read_name = bam_get_qname(aln);
    auto it = read_map.find(read_name);
    //assumption is that read pairs share a name
    if (it != read_map.end()) {
      bam1_t* mate = it->second;
      if (aln->core.flag & BAM_FREVERSE){
        merge_reads(mate, aln, amplicons);
      } else {
        merge_reads(aln, mate, amplicons);
      }
       
      //exit(0);
      //clean the mate
      bam_destroy1(mate);
      read_map.erase(it);
    } else {
      // Store the current read in the map
      read_map[read_name] = bam_dup1(aln);  // Duplicate the read to avoid overwriting
    }*/

    //TEST LINES
    //if(start_pos > 3500) continue;
    //std::string test = bam_get_qname(aln);
    //if(test != "A01535:8:HJ3YYDSX2:4:1126:24433:27305") continue;
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
      //CHANGE can do this by index
      for(uint32_t j=0; j < primers.size(); j++){
        uint32_t pstart = primers[j].get_start();
        uint32_t pend = primers[j].get_end(); 
        if (start == pstart && end == pend){          
          primers[j].add_cigarotype(cigar, aln->core.pos, nlength, seq, aux, bam_get_qname(aln), qualities, bam_is_rev(aln));
        }
      }
    }
    
  }

  std::cerr << "transforming primers" << std::endl;
  //this is super time costly
  for(uint32_t i=0; i < primers.size(); i++){
    primers[i].populate_positions();
    primers[i].transform_mutations();
  }
  std::cerr << "setting haplotypes" << std::endl;
  for (uint32_t i=0; i < primers.size(); i++){
    amplicons.set_haplotypes(primers[i]);      
  }
  /*
  //TEST LINES
  std::string amp_file = "/home/chrissy/paper_ivar_2.0/real_primer_binding/var/file_121_all.txt";
  std::ofstream file_amp(amp_file, std::ios::trunc); 
  file_amp << "POS\tALLELE\tDEPTH\tFREQ\tAMP_START\tAMP_END\n";
  file_amp.close();
  amplicons.write_out_frequencies(amp_file);
  */

  //combine amplicon counts to get total variants 
  amplicons.combine_haplotypes();
  //detect primer binding issues
  std::vector<position> variants = amplicons.variants;
  //add in primer info
  for(uint32_t i=0; i < variants.size(); i++){
    bool mutation = false;
    //establish a mutation at this position
    for(auto al : variants[i].alleles){
      if(((float)al.depth/(float)variants[i].depth) > 0.05){
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
  uint32_t test_pos = 0;
  //detect fluctuating variants
  for(uint32_t i=0; i < amplicons.max_pos; i++){
    amplicons.test_flux.clear();
    amplicons.test_test.clear();
    
    //this bit pushes all amp position vectors back to test_flux object
    amplicons.detect_abberations(i);
    if (amplicons.test_flux.size() < 2) continue;
    
    std::map<std::string, std::vector<float>> allele_maps;
    std::map<std::string, std::vector<uint32_t>> allele_depths;  
    for(uint32_t j=0; j < amplicons.test_flux.size(); j++){
      if(i == test_pos){
        std::cerr << "\n" << amplicons.test_test[j] << std::endl;
      }
      uint32_t total_depth = amplicons.test_flux[j].depth;
      //actually, we use ungapped depth
      for(auto pos : amplicons.test_flux[j].alleles){
        //if(pos.nuc == "-"){
        //  total_depth -= pos.depth;
        //}
        if(i == test_pos){
          std::cerr << pos.nuc << " " << pos.depth << std::endl;
        }
      }

      if(total_depth < 50){
         continue;
      }
      if(i == test_pos){
        std::cerr << "total depth " << total_depth << std::endl;
      }
      std::vector<allele> ad  = amplicons.test_flux[j].alleles;
      for(uint32_t k=0; k < ad.size(); k++){
        std::string nuc = ad[k].nuc;
        //if(nuc == "-") continue;
        uint32_t ad_depth = ad[k].depth;
        if(ad_depth == 0) continue;
        if(i == test_pos){
          std::cerr << nuc << " " << ad_depth << " " << total_depth << std::endl;
        }
        float t = (float)ad_depth / (float)total_depth; 
        if (allele_maps.find(nuc) != allele_maps.end()){
          if(i == test_pos){
            std::cerr << "allele map " << nuc << " " << t << std::endl;
          }
          allele_maps[nuc].push_back(t);
          allele_depths[nuc].push_back(ad_depth);
        } else {
          std::vector<float> tmp;
          std::vector<uint32_t> tmp_depths;
          tmp.push_back(t);
          tmp_depths.push_back(ad_depth);
          if(i == test_pos){
            std::cerr << "allele map " << nuc << " " << t << std::endl;
          }
          allele_maps[nuc] = tmp;
          allele_depths[nuc] = tmp_depths;
        }
      }
    }
    std::map<std::string, std::vector<float>>::iterator it;
    for (it = allele_maps.begin(); it != allele_maps.end(); it++){
      float sd = (float)calculate_standard_deviation(it->second, allele_depths[it->first]);
      if(i == test_pos){
        for (auto bit = allele_depths.begin(); bit != allele_depths.end(); ++bit) {
          std::cerr << "allele depths " << bit->first << std::endl;
          for(auto b : bit->second){
            std::cerr << b << " ";
          }
          std::cerr << "\n";
        }
        std::cerr << i << " std " << sd << " " << it->first << std::endl;
        std::cerr << "allele frequencies" << std::endl;
        for(auto x : it->second){
          std::cerr << x << std::endl;
        }
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
              }
            }
        }
        flagged_positions.push_back(i);
        break;
      }
    }
  }
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
  //file << "POS\tALLELE\tDEPTH\tFREQ\tGAPPED_FREQ\tAVG_QUAL\tFLAGGED_POS\tAMP_MASKED\tPRIMER_BINDING\n";
  file << "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\tGAPPED_FREQ\tFLAGGED_POS\tAMP_MASKED\tPRIMER_BINDING\n";
  for(uint32_t i=0; i < variants.size(); i++){
    //find the depth of the deletion to calculate upgapped depth
    float del_depth = 0;
    for(uint32_t j=0; j < variants[i].alleles.size(); j++){
      if(variants[i].alleles[j].nuc == "-"){
        del_depth += (float)variants[i].alleles[j].depth;
        break;  
      }
    }
    for(uint32_t j=0; j < variants[i].alleles.size(); j++){
      float freq = (float)variants[i].alleles[j].depth / ((float)variants[i].depth-del_depth);
      float gapped_freq = (float)variants[i].alleles[j].depth / (float)variants[i].depth;
      if((float)variants[i].alleles[j].depth == 0){
        continue;
      }
      file << "NA\t"; //region
      file << std::to_string(variants[i].pos) << "\t";
      file << "NA\t"; //ref
      std::string test_string = variants[i].alleles[j].nuc;
      file << variants[i].alleles[j].nuc << "\t";
      file << "NA\t"; //ref dp
      file << "NA\t"; //ref rv
      file << "NA\t"; //ref qual   
      file << std::to_string(variants[i].alleles[j].depth) << "\t"; //alt dp
      file << "NA\t"; //alt rv
      file << std::to_string((float)variants[i].alleles[j].mean_qual / (float)variants[i].alleles[j].depth) << "\t"; //alt qual
      file << std::to_string(freq) << "\t"; //alt freq
      file << std::to_string(variants[i].depth) << "\t"; //total dp
      file << "NA\t"; //pval
      file << "NA\t"; //pass
      file << "NA\t"; //gff feature
      file << "NA\t"; //ref codon
      file << "NA\t"; //ref aa
      file << "NA\t"; //alt codon
      file << "NA\t"; //alt aa
      file << "NA\t"; //pos aa
      file << std::to_string(gapped_freq) << "\t";
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
      file << "\n";
    } 
  } 
  file.close();
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return(retval);
}
