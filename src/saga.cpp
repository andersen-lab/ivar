#include "saga.h"
#include "ref_seq.h"
#include "parse_gff.h"
#include <fstream>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <tuple>

std::string join_uint32_vector(const std::vector<uint32_t>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec[i];
    }
    return oss.str();
}
std::string join_double_vector(const std::vector<double>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec[i];
    }
    return oss.str();
}

uint32_t calculate_reference_depth(position var, char ref){
  uint32_t depth = 0;
  for(auto allele : var.alleles){
    if(allele.nuc[0] == ref){
      depth = allele.depth;
      break;
    }
  }
  return(depth);
}
uint32_t calculate_reference_qual(position var, char ref){
  uint32_t qual= 0;
  for(auto allele : var.alleles){
    if(allele.nuc[0] == ref){
      qual = allele.mean_qual;
      break;
    }
  }
  return(qual);
}

std::vector<allele> find_deletions_next(std::vector<position> variants, uint32_t pos){
  std::vector<allele> deletions;
  for(auto var : variants){
    if(var.pos-1 == pos){
      for(uint32_t j=0; j < var.alleles.size(); j++){
        if (var.alleles[j].nuc.find("-") != std::string::npos && var.alleles[j].depth > 0){
          deletions.push_back(var.alleles[j]);
        }    
      }
    }
  }
  return(deletions);
}

void parse_cigar(const bam1_t* read1, std::vector<uint32_t> &positions, std::vector<std::string> &bases, std::vector<uint32_t> &qualities, uint32_t total_ref_pos, uint8_t min_qual, ref_antd &refantd, std::string ref_name){
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
      std::string tmp = "-";
      positions.push_back(total_ref_pos+counter);
      for(uint32_t j=0; j < len; j++){
        tmp += refantd.get_base(total_ref_pos+counter, ref_name);        
        counter++;
      }
      bases.push_back(tmp);
      qualities.push_back((uint32_t)min_qual);
    } else if(op == 1){
      std::string tmp = "+";       
      double qual_avg = 0; 
      //collect all nucs in insertions
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        tmp +=  std::string(1, nuc);
        qual_avg += qual[j];
      }
      positions.push_back(total_ref_pos+counter-1);
      bases.push_back(tmp);
      qualities.push_back((uint32_t)qual_avg/(tmp.size()-1));
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

void merge_reads(const bam1_t* read1, const bam1_t* read2, IntervalTree &amplicons, uint8_t min_qual, ref_antd &refantd, std::string ref_name){
  //pass the forward first then reverse 
  //also assumes that the forward read starts more "left" than  the reverse
  uint32_t start_reverse = read2->core.pos; 
  uint32_t start_forward = read1->core.pos;
  uint32_t end_reverse = find_sequence_end(read2);

  //record the positions and their bases
  std::vector<uint32_t> positions1;
  std::vector<std::string> bases1;
  std::vector<uint32_t> qualities1;

  std::vector<uint32_t> positions2;
  std::vector<std::string> bases2;
  std::vector<uint32_t> qualities2;
  
  //start of read, 
  parse_cigar(read1, positions1, bases1, qualities1, start_forward, min_qual, refantd, ref_name);
  parse_cigar(read2, positions2, bases2, qualities2, start_reverse, min_qual, refantd, ref_name);  

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
int preprocess_reads(std::string bam, std::string bed, std::string bam_out, std::string cmd, std::string pair_info, int32_t primer_offset, uint32_t min_depth, uint8_t min_qual, std::string ref_file){
  if(ref_file.empty()){
    std::cerr << "Please provide reference." << std::endl;
    exit(1);
  }
  //load the reference
  ref_antd refantd(ref_file, "");
  
  bool development_mode = true;
  int retval = 0;
  std::vector<primer> primers;
  if (!bed.empty()) {
    primers = populate_from_file(bed, primer_offset);
  }
  std::string gff_path = "";
  IntervalTree amplicons;
  if (!pair_info.empty() && primers.size() != 0) {
    amplicons = populate_amplicons(pair_info, primers);    
    amplicons.inOrder();
    amplicons.get_max_pos();
    amplicons.amplicon_position_pop();
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
  std::vector<bam1_t *> alns;
   

  //iterate through reads
  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  aln = bam_init1();

  uint32_t last_position=0;
  int tid=0;
  while (sam_read1(in, header, aln) >= 0) {
    uint32_t end_pos = aln->core.pos + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));
    if (end_pos > last_position) {
      last_position = end_pos;
      //get the region
      tid = aln->core.tid;
    }
  }
  const std::string ref_name = (std::string)header->target_name[tid];
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);

  in = sam_open(bam.c_str(), "r");
  header = sam_hdr_read(in);
  aln = bam_init1();
  amplicons.populate_variants(last_position);

  //hold the reads until it's mate can be found
  std::unordered_map<std::string, bam1_t*> read_map;
   
  // Iiterate through reads
  while (sam_read1(in, header, aln) >= 0) {
    //get the name of the read    
    std::string read_name = bam_get_qname(aln);

    if (!(aln->core.flag & BAM_FPAIRED) || !(aln->core.flag & BAM_FPROPER_PAIR)){
      //if the read is unpaired try to assign it to an amplicon anyways
      std::vector<uint32_t> positions;
      std::vector<std::string> bases;
      std::vector<uint32_t> qualities;
      uint32_t start_read = aln->core.pos;
      uint32_t end_read = find_sequence_end(aln);
      parse_cigar(aln, positions, bases, qualities, start_read, min_qual, refantd, ref_name);
      bool found_amplicon = false;
      uint32_t amp_dist = 429496729;
      uint32_t amp_start = 0;
      amplicons.find_read_amplicon(start_read, end_read, found_amplicon, read_name, amp_start, amp_dist);   
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
        merge_reads(mate, aln, amplicons, min_qual, refantd, ref_name);
      } else {
        merge_reads(aln, mate, amplicons, min_qual, refantd, ref_name);
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
      parse_cigar(it->second, positions, bases, qualities, start_read, min_qual, refantd, ref_name);
      bool found_amplicon = false;
      uint32_t amp_dist = 429496729;
      uint32_t amp_start = 0;
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
  uint32_t counter = 1;
  amplicons.combine_haplotypes(counter);
  //detect primer binding issues
  std::vector<position> variants = amplicons.variants;

  std::vector<uint32_t> flagged_positions; 
  std::vector<float> std_deviations;
  std::vector<std::string> pos_nuc;
  uint32_t test_pos = 0;
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
      }
    }
  }
  std::cerr << "variants size " << variants.size() << std::endl;
  std::vector<uint32_t> flagged_amplicons;
  //find amplicons that are problematic
  for(uint32_t j=0; j < variants.size(); j++){
    auto it = std::find(flagged_positions.begin(), flagged_positions.end(), variants[j].pos);
    if(it != flagged_positions.end()){
      for(auto an : variants[j].amplicon_numbers){
        auto ait = std::find(flagged_amplicons.begin(), flagged_amplicons.end(), an);
        if(ait == flagged_amplicons.end()){
          flagged_amplicons.push_back(an);
        }
      }
    }
  }
  //write variants to a file
  ofstream file;
  file.open(bam_out + ".txt", ios::trunc);
  file << "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\tGAPPED_FREQ\tGAPPED_DEPTH\tFLAGGED_POS\tAMP_MASKED\tSTD_DEV\tAMP_FREQ\tAMP_NUMBERS\n";
  for(uint32_t i=0; i < variants.size(); i++){
    if(variants[i].depth == 0) continue;
    //find the depth of the deletion to calculate upgapped depth
    double del_depth = 0;
    char ref = refantd.get_base(variants[i].pos, ref_name);
    //calculate the reference depth
    uint32_t ref_depth = calculate_reference_depth(variants[i], ref);
    double ref_qual = 0;
    if(ref_depth > 0){
      uint32_t  ref_qual_avg = calculate_reference_qual(variants[i], ref);
      ref_qual = (double)ref_qual_avg / (double)ref_depth;
    }
    uint32_t pos = variants[i].pos; 

    //deletions need to be shifted one position back
    std::vector<allele> del_alleles = find_deletions_next(variants, pos);    

    std::vector<allele> alleles = variants[i].alleles;
    //remove the deletion depth as needed
    for(uint32_t j=0; j < alleles.size(); j++){
      if(alleles[j].depth == 0) continue;
      if(alleles[j].nuc.find("-") != std::string::npos){
        del_depth += (double)alleles[j].depth;
        break;  
      }
    }
    //remove any current deletions
    std::vector<uint32_t> remove_indices;
    for(uint32_t j=0; j < alleles.size(); j++){
      if (alleles[j].nuc.find("-") != std::string::npos){
        remove_indices.push_back(j);
      } 
    }
    std::sort(remove_indices.rbegin(), remove_indices.rend());
    for(uint32_t idx : remove_indices) {
      if (idx < alleles.size()) {
        alleles.erase(alleles.begin() + idx);
      }
    }   

    //add in our bonus deletions
    if(del_alleles.size() > 0){
      alleles.insert(alleles.end(), del_alleles.begin(), del_alleles.end());
    }
    //iterate all alleles and add them in
    for(uint32_t j=0; j < alleles.size(); j++){
      if(alleles[j].depth == 0){
        continue;
      }
      /*if(alleles[j].nuc[0] == ref){
        continue;
      }*/
      double freq = (double)alleles[j].depth / ((double)variants[i].depth-del_depth);
      double gapped_freq = (double)alleles[j].depth / (double)variants[i].depth;
      file << ref_name <<"\t"; //region
      file << std::to_string(variants[i].pos) << "\t";
      file << "NA"  << "\t"; //ref
      file << alleles[j].nuc << "\t";
      file << std::to_string(ref_depth) << "\t"; //ref dp
      file << ref << "\t"; //ref rv
      file << std::to_string(ref_qual) << "\t"; //ref qual   
      file << std::to_string(alleles[j].depth) << "\t"; //alt dp
      file << "NA\t"; //alt rv
      file << std::to_string((double)alleles[j].mean_qual / (double)alleles[j].depth) << "\t"; //alt qual
      file << std::to_string(freq) << "\t"; //alt freq
      file << std::to_string(variants[i].depth-del_depth) << "\t"; //total dp
      file << "NA\t"; //pval
      file << "NA\t"; //pass
      file << "NA\t"; //gff feature
      file << "NA\t"; //ref codon
      file << "NA\t"; //ref aa
      file << "NA\t"; //alt codon
      file << "NA\t"; //alt aa
      file << "NA\t"; //pos aa
      file << std::to_string(gapped_freq) << "\t"; //gapped freq
      file << std::to_string(variants[i].depth) << "\t"; //gapped depth
      std::vector<uint32_t>::iterator it; 
      it = find(flagged_positions.begin(), flagged_positions.end(), variants[i].pos);
      if (it != flagged_positions.end()){
        file << "TRUE\t";
      } else {
        file << "FALSE\t";
      }
      bool amp_flag = false;
      for(auto an : variants[i].amplicon_numbers){
        it = find(flagged_amplicons.begin(), flagged_amplicons.end(), an);
        if(it != flagged_amplicons.end()){
          amp_flag = true;
          break;
        }
      }
      if(amp_flag) file << "TRUE\t";
      else file << "FALSE\t";
      std::string tmp = std::to_string(variants[i].pos) + alleles[j].nuc;
      std::vector<std::string>::iterator sit = find(pos_nuc.begin(), pos_nuc.end(), tmp);
      if(sit != pos_nuc.end()){
        uint32_t index = std::distance(pos_nuc.begin(), sit);
        file << std::to_string(std_deviations[index]) << "\t";
      } else {
        file << "0\t";
      }
      std::vector<double> freqs = variants[i].amplicon_frequencies[alleles[j].nuc];
      if(freqs.size() > 0){
        file << join_double_vector(freqs) << "\t";     
      } else {
        file << "NA\t";
      }
      if(variants[i].amplicon_numbers.size() > 0){
        file << join_uint32_vector(variants[i].amplicon_numbers);
      } else {
        file << "NA";
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
