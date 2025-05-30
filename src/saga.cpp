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

void calculate_reference_depth(genomic_position var, char ref, uint32_t &depth){
  for(auto allele : var.alleles){
    if(allele.nuc[0] == ref){
      depth = allele.depth;
      break;
    }
  }
}

void calculate_reference_qual(genomic_position var, char ref, uint32_t &qual){
  for(auto allele : var.alleles){
    if(allele.nuc[0] == ref){
      qual = allele.mean_qual;
      break;
    }
  }
}

std::vector<allele> find_deletions_next(genomic_position position){
  std::vector<allele> deletions;
  if(position.depth == position.gapped_depth) return(deletions);
  for(uint32_t j=0; j < position.alleles.size(); j++){
    if (position.alleles[j].nuc.find("-") != std::string::npos){
      deletions.push_back(position.alleles[j]);
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
  uint32_t mqual = (uint32_t) min_qual;

  static const char seq_nt_lookup[16] = {
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
  };


  for (uint32_t i = 0; i < read1->core.n_cigar; i++){
    uint32_t op = bam_cigar_op(cigar1[i]);
    uint32_t len = bam_cigar_oplen(cigar1[i]);
    uint32_t counter = 0;
    if(op == 0){
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        uint32_t tqual = qual[j];
        if(tqual < mqual){
          counter++;
          continue;
        }
        //char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        uint8_t base_code = bam_seqi(seq_field1, j);
        char nuc = seq_nt_lookup[base_code];

        std::string tmp = std::string(1, nuc);
        positions.push_back(total_ref_pos+counter);
        bases.push_back(tmp);
        qualities.push_back(qual[j]);
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
      qualities.push_back(min_qual);
    } else if(op == 1){
      std::string tmp = "+";
      double qual_avg = 0;
      //collect all nucs in insertions
      for(uint32_t j=total_query_pos; j < total_query_pos+len; j++){
        //char nuc = seq_nt16_str[bam_seqi(seq_field1, j)];
        uint8_t base_code = bam_seqi(seq_field1, j);
        char nuc = seq_nt_lookup[base_code];
        tmp +=  std::string(1, nuc);
        qual_avg += qual[j];
      }
      if(((uint32_t)qual_avg/(tmp.size()-1)) < min_qual) continue;
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

void merge_reads(const bam1_t* read1, const bam1_t* read2, IntervalTree &amplicons, uint8_t min_qual, ref_antd &refantd, std::string ref_name, std::vector<genomic_position> &global_positions){
  //pass the forward first then reverse
  //also assumes that the forward read starts more "left" than  the reverse
  uint32_t start_reverse = read2->core.pos;
  uint32_t start_forward = read1->core.pos;
  uint32_t end_reverse = bam_endpos(read2);
  //record the positions and their bases
  std::vector<uint32_t> positions1, positions2;
  std::vector<std::string> bases1, bases2;
  std::vector<uint32_t> qualities1, qualities2;

  parse_cigar(read1, positions1, bases1, qualities1, start_forward, min_qual, refantd, ref_name);
  parse_cigar(read2, positions2, bases2, qualities2, start_reverse, min_qual, refantd, ref_name);

  // Index positions to avoid repeated scanning
  std::unordered_map<uint32_t, std::vector<uint32_t>> pos1_indices, pos2_indices;
  for (uint32_t i = 0; i < positions1.size(); ++i)
    pos1_indices[positions1[i]].push_back(i);
  for (uint32_t i = 0; i < positions2.size(); ++i)
    pos2_indices[positions2[i]].push_back(i);

  //find all unique positions we need to cover
  std::unordered_set<uint32_t> unique_elements(positions1.begin(), positions1.end());
  unique_elements.insert(positions2.begin(), positions2.end());
  std::vector<uint32_t> final_positions;
  std::vector<uint32_t> final_qualities;
  std::vector<std::string> final_bases;

  //for every position make sure the bases and qualities match
  //make sure insertions match
  for(auto pos : unique_elements){
    auto it1 = pos1_indices.find(pos);
    auto it2 = pos2_indices.find(pos);

    bool in1 = (it1 != pos1_indices.end());
    bool in2 = (it2 != pos2_indices.end());

    if (in1 && !in2) {
      for (uint32_t idx : it1->second) {
        final_positions.push_back(positions1[idx]);
        final_qualities.push_back(qualities1[idx]);
        final_bases.push_back(bases1[idx]);
      }
    } else if (!in1 && in2) {
      for (uint32_t idx : it2->second) {
        final_positions.push_back(positions2[idx]);
        final_qualities.push_back(qualities2[idx]);
        final_bases.push_back(bases2[idx]);
      }
    } else if (in1 && in2) {
      const auto& idxs1 = it1->second;
      const auto& idxs2 = it2->second;

      if (idxs1.size() != idxs2.size()) continue; // Skip if mismatch in count

      bool match = true;
      for (size_t j = 0; j < idxs1.size(); ++j) {
        if (bases1[idxs1[j]] != bases2[idxs2[j]]) {
          match = false;
          break;
        }
      }

      if (match) {
        for (size_t j = 0; j < idxs1.size(); ++j) {
          final_positions.push_back(positions1[idxs1[j]]);
          final_qualities.push_back(qualities1[idxs1[j]]);
          final_bases.push_back(bases1[idxs1[j]]);
        }
      }
    }
  }

  //find assigned amplicon and populate position vector
  uint32_t amp_dist = 429496729;
  ITNode *node=NULL;
  amplicons.find_read_amplicon(start_forward, end_reverse, node, amp_dist);
  if(node == NULL){
    add_variants(final_positions, final_bases, final_qualities, global_positions);
  } else {
    assign_read(node, final_positions, final_bases, final_qualities, global_positions);
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

double  calculate_standard_deviation_weighted(std::vector<double> values, std::vector<uint32_t> weights) {
    double weighted_sum = 0.0, total_weight = 0.0;

    // Compute weighted mean
    for (size_t i = 0; i < values.size(); ++i) {
        weighted_sum += values[i] * weights[i];
        total_weight += weights[i];
    }
    double mean = weighted_sum / total_weight;

    // Compute weighted variance
    double variance = 0.0f;
    for (size_t i = 0; i < values.size(); ++i) {
        variance += weights[i] * std::pow(values[i] - mean, 2);
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

  //initate variants
  std::vector<genomic_position> global_positions;
  //populate blank vector
  if(amplicons.max_pos != 0){
    populate_positions(global_positions, amplicons.max_pos);
  }
  //assign overlap regions
  amplicons.calculate_overlaps(global_positions);

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
      uint32_t end_read = bam_endpos(aln);
      parse_cigar(aln, positions, bases, qualities, start_read, min_qual, refantd, ref_name);
      uint32_t amp_dist = 429496729;
      ITNode *node=NULL;
      amplicons.find_read_amplicon(start_read, end_read, node, amp_dist);
      if(node == NULL){
        add_variants(positions, bases, qualities, global_positions);
      } else{
        assign_read(node, positions, bases, qualities, global_positions);
      }
      continue;
    }
    auto it = read_map.find(read_name);
    //assumption is that read pairs share a name
    //execute if we've already seen the mate
    if (it != read_map.end()) {
      bam1_t* mate = it->second;
      if (aln->core.flag & BAM_FREVERSE){
        merge_reads(mate, aln, amplicons, min_qual, refantd, ref_name, global_positions);
      } else {
        merge_reads(aln, mate, amplicons, min_qual, refantd, ref_name, global_positions);
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
      std::vector<uint32_t> positions;
      std::vector<std::string> bases;
      std::vector<uint32_t> qualities;
      uint32_t start_read = it->second->core.pos;
      uint32_t end_read = bam_endpos(it->second);
      parse_cigar(it->second, positions, bases, qualities, start_read, min_qual, refantd, ref_name);
      uint32_t amp_dist = 429496729;
      ITNode *node=NULL;
      amplicons.find_read_amplicon(start_read, end_read, node, amp_dist);
      if(node == NULL){
        add_variants(positions, bases, qualities, global_positions);
      } else{
        assign_read(node, positions, bases, qualities, global_positions);
      }
  }

  //combine amplicon counts to get total variants
  uint32_t counter = 1;
  combine_haplotypes(global_positions);
  //TODO build in quality filters, and min depth filter
  std::vector<ITNode*> flagged_amplicons = calculate_amplicon_variation(global_positions, min_depth, min_qual);
  set_amplicon_flag(flagged_amplicons, global_positions);
  //write variants to a file
  std::unordered_map<std::string, std::vector<double>> allele_frequencies;
  std::vector<allele> del_alleles;
  std::vector<uint32_t> amplicon_numbers;
  std::string blank = "NA\t";

  ofstream file;
  file.open(bam_out + ".txt", ios::trunc);
  file << "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\tGAPPED_FREQ\tGAPPED_DEPTH\tFLAGGED_POS\tAMP_MASKED\tSTD_DEV\tAMP_FREQ\tAMP_NUMBERS\n";
  for(auto var : global_positions){
    if(var.depth == 0) continue;

    char ref = refantd.get_base(var.pos, ref_name);
    //calculate the reference depth
    uint32_t ref_depth=0, ref_qual=0, ref_qual_avg=0;
    calculate_reference_depth(var, ref, ref_depth);
    if(ref_depth > 0){
      calculate_reference_qual(var, ref, ref_qual_avg);
      ref_qual = (double)ref_qual_avg / (double)ref_depth;
    }
    //deletions need to be shifted one position back
    del_alleles.clear();
    if(var.pos+1 < global_positions.size()){
      del_alleles = find_deletions_next(global_positions[var.pos+1]);
    }

    //remove any current deletions
    var.alleles.erase(
    std::remove_if(var.alleles.begin(), var.alleles.end(), [](const allele &a) {
      return a.nuc.find('-') != std::string::npos;
    }),var.alleles.end());

    //add in our bonus deletions
    if(del_alleles.size() > 0){
      var.alleles.insert(var.alleles.end(), del_alleles.begin(), del_alleles.end());
    }
    //get amplicon specific frequencies
    allele_frequencies.clear();
    collect_allele_frequencies(var.amplicons, allele_frequencies);
    //get the amplicon numbers
    amplicon_numbers.clear();
    get_amplicon_numbers(var.amplicons, amplicon_numbers);

    //iterate all alleles and add them in
    for(uint32_t j=0; j < var.alleles.size(); j++){
      if(var.alleles[j].depth == 0){
        continue;
      }
      double freq = (double)var.alleles[j].depth / ((double)var.depth);
      double gapped_freq = (double)var.alleles[j].depth / (double)var.gapped_depth;
      file << ref_name <<"\t"; //region
      file << std::to_string(var.pos) << "\t";
      file << blank; //ref
      file << var.alleles[j].nuc << "\t";
      file << std::to_string(ref_depth) << "\t"; //ref dp
      file << ref << "\t"; //ref rv
      file << std::to_string(ref_qual) << "\t"; //ref qual
      file << std::to_string(var.alleles[j].depth) << "\t"; //alt dp
      file << blank; //alt rv
      file << std::to_string((double)var.alleles[j].mean_qual / (double)var.alleles[j].depth) << "\t"; //alt qual
      file << std::to_string(freq) << "\t"; //alt freq
      file << std::to_string(var.depth) << "\t"; //total dp ungapped
      file << blank; //pval
      file << blank; //pass
      file << blank; //gff feature
      file << blank; //ref codon
      file << blank; //ref aa
      file << blank; //alt codon
      file << blank; //alt aa
      file << blank; //pos aa
      file << std::to_string(gapped_freq) << "\t"; //gapped freq
      file << std::to_string(var.gapped_depth) << "\t"; //gapped depth
      //handle variant level fluctuation
      if (var.flux){
        file << "TRUE\t";
      } else {
        file << "FALSE\t";
      }
      //handle amplicon level fluctuation
      if(var.amp_flux) file << "TRUE\t";
      else file << "FALSE\t";
      //standard deviation across amplicons
      file << std::to_string(var.alleles[j].stddev) << "\t";
      //the amplicon-specific frequencies of the allele
      std::vector<double> freqs = allele_frequencies[var.alleles[j].nuc];
      if(freqs.size() > 0){
        file << join_double_vector(freqs) << "\t";
      } else {
        file << blank;
      }

      //write out amplicons it's assigned to
      if(amplicon_numbers.size() > 0){
        file << join_uint32_vector(amplicon_numbers);
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
