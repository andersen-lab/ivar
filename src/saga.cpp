#include "saga.h"
#include <fstream>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <tuple>

void parse_cigar(const bam1_t* read1, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities, uint32_t total_ref_pos, uint8_t min_qual){
  uint32_t total_query_pos=0;
  const uint8_t* seq_field1 = bam_get_seq(read1);
  uint32_t *cigar1 = bam_get_cigar(read1);
  uint8_t* qual = bam_get_qual(read1);
  total_ref_pos += 1;
  for (uint32_t i = 0; i < read1->core.n_cigar; i++){
    uint32_t op = bam_cigar_op(cigar1[i]);
    uint32_t len = bam_cigar_oplen(cigar1[i]);    
    uint32_t counter = 0;
    if(op == 0){ 
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        std::string tmp = std::string(1, nuc);
        positions.push_back(total_ref_pos+counter);
        bases.push_back(tmp);
        qualities.push_back((uint32_t)qual[j]);
        counter++;
      }
    } else if(op == 2){    
      for(uint32_t j=0; j < len; j++){
        positions.push_back(total_ref_pos+counter);
        bases.push_back("-");
        qualities.push_back((uint32_t)min_qual);
        counter++;
      }
    } else if(op == 1){
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        std::string tmp = "+" + std::string(1, nuc);
        positions.push_back(total_ref_pos+counter);
        bases.push_back(tmp);
        qualities.push_back((uint32_t)qual[j]);
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

void merge_reads(const bam1_t* read1, const bam1_t* read2, IntervalTree &amplicons, uint8_t min_qual){
  //pass the forward first then reverse 
  //also assumes that the forward read starts more "left" than  the reverse
  //std::cerr << "merging reads"<< std::endl;  
  uint32_t start_reverse = read2->core.pos; 
  uint32_t start_forward = read1->core.pos;
  uint32_t end_reverse = find_sequence_end(read2);

  //std::cerr << start_forward << "," << end_reverse << std::endl;
  
  //record the positions and their bases
  std::vector<uint32_t> positions1;
  std::vector<std::string> bases1;
  std::vector<uint32_t> qualities1;

  std::vector<uint32_t> positions2;
  std::vector<std::string> bases2;
  std::vector<uint32_t> qualities2;
  
  //start of read, 
  parse_cigar(read1, positions1, bases1, qualities1, start_forward, min_qual);
  parse_cigar(read2, positions2, bases2, qualities2, start_reverse, min_qual);  

  //find all unique positions we need to cover
  std::unordered_set<uint32_t> unique_elements(positions1.begin(), positions1.end());
  unique_elements.insert(positions2.begin(), positions2.end());

  std::vector<uint32_t> final_positions;
  std::vector<uint32_t> final_qualities;
  std::vector<std::string> final_bases;

  //for every position make sure the bases and qualities match
  //make sure insertions match
  for(auto pos : unique_elements){
    auto first_it = std::find(positions1.begin(), positions1.end(), pos);
    auto second_it = std::find(positions2.begin(), positions2.end(), pos);
    
    //if its in the first read but not the second
    if(first_it != positions1.end() && second_it == positions2.end()){
      auto it = positions1.begin();
      while ((it = std::find(it, positions1.end(), pos)) != positions1.end()) {
        uint32_t index1 = std::distance(positions1.begin(), it);
        final_positions.push_back(positions1[index1]);
        final_qualities.push_back(qualities1[index1]);
        final_bases.push_back(bases1[index1]);
        ++it;
      }
    } else if(first_it == positions1.end() && second_it != positions2.end()){
      //if its in the second read but not the first
      auto it = positions2.begin();
      while ((it = std::find(it, positions2.end(), pos)) != positions2.end()) {
        uint32_t index2 = std::distance(positions2.begin(), it);
        final_positions.push_back(positions2[index2]);
        final_qualities.push_back(qualities2[index2]);
        final_bases.push_back(bases2[index2]);
        ++it;
      }
    } else if(first_it != positions1.end() && second_it != positions2.end()){
      //if its in both reads
      //check if the number of times the position occurs is the same
      uint32_t count1 = std::count(positions1.begin(), positions1.end(), pos);
      uint32_t count2 = std::count(positions2.begin(), positions2.end(), pos);
      if(count1 != count2){continue;}
   
      std::vector<uint32_t> tmp_pos2; 
      std::vector<std::string> tmp_base2; 
      std::vector<uint32_t> tmp_qual2; 
      //go get all the second read info
      auto sit = positions2.begin();
      while ((sit = std::find(sit, positions2.end(), pos)) != positions2.end()) {
        uint32_t index2 = std::distance(positions2.begin(), sit);
        tmp_pos2.push_back(positions2[index2]);
        tmp_qual2.push_back(qualities2[index2]);
        tmp_base2.push_back(bases2[index2]);
        ++sit;
      }
      std::vector<uint32_t> tmp_pos1; 
      std::vector<std::string> tmp_base1; 
      std::vector<uint32_t> tmp_qual1; 
      //go get all the first read info
      auto fit = positions1.begin();
      while ((fit = std::find(fit, positions1.end(), pos)) != positions1.end()) {
        uint32_t index1 = std::distance(positions1.begin(), fit);
        tmp_pos1.push_back(positions1[index1]);
        tmp_qual1.push_back(qualities1[index1]);
        tmp_base1.push_back(bases1[index1]);
        ++fit;
      }
      //now do the base and quality comparison
      bool use= true;
      for(uint32_t j=0; j < tmp_pos1.size(); j++){
        if(tmp_base1[j] != tmp_base2[j]){
          use = false;
          break;
        }
      }
      if(use) {
        for(uint32_t j=0; j < tmp_pos1.size(); j++){
          final_positions.push_back(tmp_pos1[j]);
          final_qualities.push_back(tmp_qual1[j]);
          final_bases.push_back(tmp_base1[j]);
        }
      }
    }
  }

  //find assigned amplicon and populate position vector
  bool found_amplicon = false;
  uint32_t amp_dist = 429496729;
  uint32_t amp_start = 0;
  //std::cerr << bam_get_qname(read1) << std::endl;
  amplicons.find_read_amplicon(start_forward, end_reverse, found_amplicon, bam_get_qname(read1), amp_start, amp_dist);   
  if(!found_amplicon){
    amplicons.add_read_variants(final_positions, final_bases, final_qualities, min_qual);
  } else {
    amplicons.assign_read_amplicon(amp_start, final_positions, final_bases, final_qualities, min_qual);
  }
}

void generate_range(uint32_t start, uint32_t end, std::vector<uint32_t> &result) {
  if (start > end) {
    std::swap(start, end);
  }
  for (uint32_t i = start; i <= end; ++i) {
    result.push_back(i);
  }
}

double calculate_standard_deviation(std::vector<double> data) {
    double sum = 0, mean;
    mean = std::accumulate(data.begin(), data.end(), 0.0f) / data.size();
    for (double val : data) {
        sum += std::pow(val - mean, 2);
    }
    return std::sqrt(sum / data.size());
}

float calculate_standard_deviation_weighted(std::vector<float> values, std::vector<uint32_t> weights) {
    float weighted_sum = 0.0f, total_weight = 0.0f;

    // Compute weighted mean
    for (size_t i = 0; i < values.size(); ++i) {
        weighted_sum += values[i] * (float)weights[i];
        total_weight += (float)weights[i];
    }
    float mean = weighted_sum / total_weight;

    // Compute weighted variance
    float variance = 0.0f;
    for (size_t i = 0; i < values.size(); ++i) {
        variance += (float)weights[i] * std::pow(values[i] - mean, 2);
    }
    variance /= total_weight;

    return std::sqrt(variance);
}

//first main function call
int preprocess_reads(std::string bam, std::string bed, std::string bam_out, std::string cmd, std::string pair_info, int32_t primer_offset, uint32_t min_depth, uint8_t min_qual){
  bool calculate_amplicons = true; 
  bool development_mode = true;
  int retval = 0;
  std::vector<primer> primers;
  if (!bed.empty()) {
    primers = populate_from_file(bed, primer_offset);
    if (primers.size() == 0) {
      calculate_amplicons = false;
      //std::cerr << "Exiting." << std::endl;
      //return -1;
    }
  }
  std::string gff_path = "";

  //int max_primer_len = get_bigger_primer(primers);
  // get coordinates of each amplicon
  IntervalTree amplicons;
  if (!pair_info.empty() && calculate_amplicons) {
    amplicons = populate_amplicons(pair_info, primers);    
    amplicons.inOrder();
    amplicons.get_max_pos();
    amplicons.amplicon_position_pop();
    std::cerr << "Maximum position " << amplicons.max_pos << std::endl;
    if(amplicons.max_pos == 0) calculate_amplicons = false;
  } else{
    //std::cerr << "Exiting." << std::endl;
    //return -1;
    calculate_amplicons = false;
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


  // Iterate through reads
  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  aln = bam_init1();

  uint32_t last_position=0;
  while (sam_read1(in, header, aln) >= 0) {
    uint32_t end_pos = aln->core.pos + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));
    if (end_pos > last_position) {
      last_position = end_pos;
    }
  }
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);

  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  aln = bam_init1();
  amplicons.populate_variants(last_position);

  //hold the reads until it's mate can be found
  std::unordered_map<std::string, bam1_t*> read_map;
   
  // Iiterate through reads
  while (sam_read1(in, header, aln) >= 0) {
    //get the name of the read    
    char strand = '+';
    uint32_t start_pos = -1;
    std::string read_name = bam_get_qname(aln);
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
    uint8_t *qualities = bam_get_qual(r);
    if(!calculate_amplicons){
      //TODO
      //amplicons.add_read_variants(cigar, aln->core.pos, nlength, seq, aux, qualities, bam_get_qname(aln), min_qual);
      continue;
    }
    if (!(aln->core.flag & BAM_FPAIRED) || !(aln->core.flag & BAM_FPROPER_PAIR)){
      //if the read is unpaired try to assign it to an amplicon anyways
      std::vector<uint32_t> positions;
      std::vector<std::string> bases;
      std::vector<uint32_t> qualities;
      uint32_t start_read = aln->core.pos;
      uint32_t end_read = find_sequence_end(aln);
      parse_cigar(aln, positions, bases, qualities, start_read, min_qual);
      bool found_amplicon = false;
      uint32_t amp_dist = 429496729;
      uint32_t amp_start = 0;
      amplicons.find_read_amplicon(start_read, end_read, found_amplicon, read_name, amp_start, amp_dist);   
      //std::cerr << bam_get_qname(aln) << std::endl;
      if(!found_amplicon){
        amplicons.add_read_variants(positions, bases, qualities, min_qual);
      } else{
        amplicons.assign_read_amplicon(amp_start, positions, bases, qualities, min_qual);
      }
      continue;
    }
    auto it = read_map.find(read_name);
    //assumption is that read pairs share a name
    //execute if we've already seen the mate
    if (it != read_map.end()) {
      bam1_t* mate = it->second;
      if (aln->core.flag & BAM_FREVERSE){
        merge_reads(mate, aln, amplicons, min_qual);
      } else {
        merge_reads(aln, mate, amplicons, min_qual);
      }
      //clean the mate
      bam_destroy1(mate);
      read_map.erase(it);
    } else {
      // Store the current read in the map
      read_map[read_name] = bam_dup1(aln);  // Duplicate the read to avoid overwriting
    }
  }
  //for reads that aren't flagged as unmapped but are for some reason
    for(auto it = read_map.begin(); it != read_map.end(); ++it){
      std::string read_name = bam_get_qname(it->second);    
      std::vector<uint32_t> positions;
      std::vector<std::string> bases;
      std::vector<uint32_t> qualities;
      uint32_t start_read = it->second->core.pos;
      uint32_t end_read = find_sequence_end(it->second);
      parse_cigar(it->second, positions, bases, qualities, start_read, min_qual);
      bool found_amplicon = false;
      uint32_t amp_dist = 429496729;
      uint32_t amp_start = 0;
      //std::cerr << read_name << std::endl;
      amplicons.find_read_amplicon(start_read, end_read, found_amplicon, read_name, amp_start, amp_dist); 
      if(!found_amplicon){
        amplicons.add_read_variants(positions, bases, qualities, min_qual);
      } else{
        amplicons.assign_read_amplicon(amp_start, positions, bases, qualities, min_qual);
      }
  }

  //TEST LINES
  if(development_mode){
    std::string amp_file = bam_out + "_all.txt";
    std::ofstream file_amp(amp_file, std::ios::trunc); 
    file_amp << "POS\tALLELE\tDEPTH\tFREQ\tAMP_START\tAMP_END\n";
    file_amp.close();
    amplicons.write_out_frequencies(amp_file);
  }

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
  std::vector<uint32_t> flagged_positions; 
  std::vector<float> std_deviations;
  std::vector<std::string> pos_nuc;
  std::vector<std::string> freq_strings;
  uint32_t test_pos = 8007;
  //detect fluctuating variants across amplicons
  for(uint32_t i=0; i < amplicons.max_pos; i++){
    amplicons.test_flux.clear();
    amplicons.test_test.clear();
    
    //this bit pushes all amp position vectors back to test_flux object
    amplicons.detect_abberations(i);
    if (amplicons.test_flux.size() < 2) continue;
    
    std::map<std::string, std::vector<float>> allele_maps;
    std::vector<uint32_t> amplicon_depths;
    for(uint32_t j=0; j < amplicons.test_flux.size(); j++){
      if(i == test_pos){
        std::cerr << "\n" << amplicons.test_test[j] << std::endl;
      }
      //this is the total depth for the amplicon!
      uint32_t total_depth = amplicons.test_flux[j].depth;
      for(auto pos : amplicons.test_flux[j].alleles){
        //remove things with crappy quality
        if(pos.mean_qual < min_qual && pos.nuc != "-"){
          total_depth -= pos.depth;
        }
        if(i == test_pos){
          std::cerr << pos.nuc << " " << pos.depth << std::endl;
        }
      }

      if(total_depth < min_depth){
         continue;
      }
      if(i == test_pos){
        std::cerr << "total depth " << total_depth << std::endl;
      }
      amplicon_depths.push_back(total_depth);
      std::vector<allele> ad  = amplicons.test_flux[j].alleles;
      for(uint32_t k=0; k < ad.size(); k++){
        if(ad[k].mean_qual < min_qual && ad[k].nuc != "-") continue;
        std::string nuc = ad[k].nuc;
        uint32_t ad_depth = ad[k].depth;
        if(ad_depth == 0) continue;
        if(i == test_pos){
          std::cerr << nuc << " " << ad_depth << " " << total_depth << std::endl;
        }
        float t = (float)ad_depth / (float)total_depth; 
        if (allele_maps.find(nuc) != allele_maps.end()){
          allele_maps[nuc].push_back(t);
        } else {
          std::vector<float> tmp;
          std::vector<uint32_t> tmp_depths;
          tmp.push_back(t);
          tmp_depths.push_back(ad_depth);
          allele_maps[nuc] = tmp;
        }
      }
    }
    std::map<std::string, std::vector<float>>::iterator it;
    for (it = allele_maps.begin(); it != allele_maps.end(); it++){
      if(it->second.size() == 1) continue; //we don't standard dev one thing
      float sd = (float)calculate_standard_deviation_weighted(it->second, amplicon_depths);
      //float sd = (float)calculate_standard_deviation(it->second);
      if(i == test_pos){
        std::cerr << i << " std " << sd << " " << it->first << std::endl;
        std::cerr << "allele frequencies" << std::endl;
        for(auto x : it->second){
          std::cerr << x << " ";
        }
        std::cerr << "\n";
      }
      //TODO this is hard coded, consider it
      if (sd >= 0.03){
        flagged_positions.push_back(i);
        std::string tmp = std::to_string(i) + it->first;
        pos_nuc.push_back(tmp);
        std_deviations.push_back(sd);
        //save the amplicon specific frequencies
        std::ostringstream oss;
        for (uint32_t k = 0; k < it->second.size(); ++k) {
          oss << it->second[k];
          if (k != it->second.size() - 1) {
            oss << ",";
          }
        }
        freq_strings.push_back(oss.str());
        //break;
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


  //assign positions to amplicons
  //first find the biggest pos covered by variants
  uint32_t max_pos = variants.size();
  for(uint32_t i=0; i < max_pos; i++){   
      std::vector<uint32_t> overlaps;
      uint32_t counter = 1;
      amplicons.detect_position_amplicons(i, counter, overlaps);
      std::stringstream ss;
      for (uint32_t j = 0; j < overlaps.size(); ++j) {
        ss << overlaps[j];
        if (j < overlaps.size() - 1) {
          ss << ",";
        }
      }
      if(ss.str() != ""){
        variants[i].amplicon_numbers = ss.str();
      }
  }

  //write variants to a file
  ofstream file;
  file.open(bam_out + ".txt", ios::trunc);
  //file << "POS\tALLELE\tDEPTH\tFREQ\tGAPPED_FREQ\tAVG_QUAL\tFLAGGED_POS\tAMP_MASKED\tPRIMER_BINDING\n";
  file << "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\tGAPPED_FREQ\tFLAGGED_POS\tAMP_MASKED\tPRIMER_BINDING\tSTD_DEV\tAMP_FREQ\tAMP_NUMBERS\n";
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
      std::string tmp = std::to_string(variants[i].pos) + variants[i].alleles[j].nuc;
      std::vector<std::string>::iterator sit = find(pos_nuc.begin(), pos_nuc.end(), tmp);
      if(sit != pos_nuc.end()){
        uint32_t index = std::distance(pos_nuc.begin(), sit);
        file << std::to_string(std_deviations[index]) << "\t";
        file << freq_strings[index] << "\t"; 
      } else {
        file << "0\t";
        file << "0\t";
      }
      file << variants[i].amplicon_numbers;
      file << "\n";
    } 
  } 
  file.close();
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return(retval);
}
