#include "trim_primer_quality.h"

#define round_int(x,total) ((int) (0.5 + ((float)x / float(total)) * 10000))/(float)100

int32_t get_pos_on_query(uint32_t *cigar, uint32_t ncigar, int32_t pos, int32_t ref_start){
  int cig;
  int32_t n;
  int32_t ql = 0, rl = ref_start;
  for (uint32_t i = 0; i < ncigar; ++i){
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if (bam_cigar_type(cig) & 2) { // Reference consuming
      if (pos <= rl + n) {
	if (bam_cigar_type(cig) & 1) // Query consuming
	  ql += (pos - rl);	   // n consumed reference, check if it consumes query too.
	return ql;
      }
      rl += n;
    }
    if (bam_cigar_type(cig) & 1) // Query consuming
      ql += n;
  }
  return ql;
}

// Number of bases from 3' end for reverse reads
int32_t get_pos_on_reference(uint32_t *cigar, uint32_t ncigar, uint32_t pos, uint32_t ref_start){
  int cig;
  int32_t n;
  uint32_t ql = 0, rl = ref_start;
  for (uint32_t i = 0; i < ncigar; ++i){
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if (bam_cigar_type(cig) & 1) { // Only query consuming
      if (pos <= ql + n) {
	if (bam_cigar_type(cig) & 2) // Only reference consuming
	  rl += (pos - ql);	   // n consumed reference, check if it consumes query too.
	return rl;
      }
      ql += n;
    }
    if (bam_cigar_type(cig) & 2) // Only reference consuming
      rl += n;
  }
  return rl;
}

void reverse_qual(uint8_t *q, int l){
  for (int i = 0; i < l/2; ++i){
    q[i]^=q[l-i-1];
    q[l-i-1]^=q[i];
    q[i]^=q[l-i-1];
  }
}

void reverse_cigar(uint32_t *cigar, int l){
  for (int i = 0; i < l/2; ++i){
    cigar[i]^=cigar[l-i-1];
    cigar[l-i-1]^=cigar[i];
    cigar[i]^=cigar[l-i-1];
  }
}

double mean_quality(uint8_t *a, int s, int e){
  double m = 0;
  for (int i = s; i < e; ++i){
    m += (double)a[i];
  }
  m = m/(e-s);
  return m;
}

cigar_ quality_trim(bam1_t* r, uint8_t qual_threshold, uint8_t sliding_window){
  bool reverse = false;
  uint32_t *ncigar = (uint32_t*) malloc(sizeof(uint32_t) * (r->core.n_cigar + 1)), // Maximum edit is one more element with soft mask
    *cigar = bam_get_cigar(r);
  uint8_t *qual = bam_get_qual(r);
  int32_t start_pos;
  if(((r->core.flag&BAM_FPAIRED) != 0) && bam_is_rev(r)){
    reverse = true;
    reverse_qual(qual, r->core.l_qseq);
  }
  double m = 60;
  int del_len, cig, temp;
  uint32_t i = 0, j = 0;
  cigar_ t;
  init_cigar(&t);
  if(0 > r->core.l_qseq - sliding_window)
    sliding_window = (uint32_t)r->core.l_qseq;
  while(i < (uint32_t)r->core.l_qseq){
    m = mean_quality(qual, i, i+sliding_window);
    if(m < qual_threshold)
      break;
    i++;
    if(i > (uint32_t)r->core.l_qseq - sliding_window)
      sliding_window--;
  }
  // Reverse qual back.
  if(reverse){
    reverse_qual(qual, r->core.l_qseq);
  }
  del_len = r->core.l_qseq - i;
  start_pos = get_pos_on_reference(cigar, r->core.n_cigar, del_len, r->core.pos); // For reverse reads need to set core->pos.
  if(reverse && start_pos <= r->core.pos) {
    free(ncigar);
    t.cigar = cigar;
    t.free_cig = false;
    t.nlength = r->core.n_cigar;
    t.start_pos = r->core.pos;
    return t;
  }
  int32_t n;
  i = 0;
  if(reverse){
    reverse_cigar(cigar, r->core.n_cigar);
  }
  reverse_cigar(cigar, r->core.n_cigar); // Reverse cigar and trim the beginning of read.
  while(i < r->core.n_cigar){
    if (del_len == 0){
      ncigar[j] = cigar[i];
      i++;
      j++;
      continue;
    }
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if ((bam_cigar_type(cig) & 1)){ // Consumes Query
      if(del_len >= n ){
	ncigar[j] = bam_cigar_gen(n, BAM_CSOFT_CLIP);
      } else if (del_len < n){
	ncigar[j] = bam_cigar_gen(del_len, BAM_CSOFT_CLIP);
      }
      j++;
      temp = n;
      n = std::max(n - del_len, 0);
      del_len = std::max(del_len - temp, 0);
      if(n > 0){
	ncigar[j] = bam_cigar_gen(n, cig);
	j++;
      }
    }
    i++;
  }
  reverse_cigar(ncigar, j);	// Reverse Back
  if(reverse){
    reverse_cigar(ncigar, j);
  }
  t.cigar = ncigar;
  t.nlength = j;
  t.free_cig = true;
  t.start_pos = start_pos;
  return t;
}

void print_cigar(uint32_t *cigar, int nlength){
  for (int i = 0; i < nlength; ++i){
    std::cout << ((cigar[i]) & BAM_CIGAR_MASK);
    std::cout << "-" << ((cigar[i]) >> BAM_CIGAR_SHIFT) << " ";
  }
  std::cout << std::endl;
}

cigar_ primer_trim(bam1_t *r, bool &isize_flag, int32_t new_pos, bool unpaired_rev = false){
  uint32_t *ncigar = (uint32_t*) malloc(sizeof(uint32_t) * (r->core.n_cigar + 1)), // Maximum edit is one more element with soft mask
    *cigar = bam_get_cigar(r);
  uint32_t i = 0, j = 0;
  int max_del_len = 0, cig, temp, del_len = 0;
  bool reverse = false;
  if((r->core.flag&BAM_FPAIRED) != 0 && isize_flag){ // If paired and isize > read length
    if (bam_is_rev(r)){ // If -ve strand (?)
      max_del_len = bam_cigar2qlen(r->core.n_cigar, bam_get_cigar(r)) - get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos) - 1;
      reverse_cigar(cigar, r->core.n_cigar);
      reverse = true;
    } else {
      max_del_len = get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos);
    }
  } else {			// trim without considering pairing
    if(unpaired_rev){
      max_del_len = bam_cigar2qlen(r->core.n_cigar, bam_get_cigar(r)) - get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos) - 1;
      reverse_cigar(cigar, r->core.n_cigar);
      reverse = true;
    } else {
      max_del_len = get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos);
    }
  }
  max_del_len = (max_del_len > 0) ? max_del_len : 0; // For cases where reads spans only primer region
  int32_t n, start_pos = 0, ref_add = 0;
  bool pos_start = false;
  del_len = max_del_len;
  while(i < r->core.n_cigar){
    if (del_len == 0 && pos_start){ // No more bases on query to soft clip
      ncigar[j] = cigar[i];
      i++;
      j++;
      continue;
    }
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if(del_len ==0 && (bam_cigar_type(cig) & 1) && (bam_cigar_type(cig) & 2)){ // After soft clipping of query complete, keep incrementing start_pos until first base that consumes both query and ref
      pos_start = true;
      continue;
    }
    ref_add = n;
    if ((bam_cigar_type(cig) & 1)){ // Consumes Query
      if(del_len >= n ){
	ncigar[j] = bam_cigar_gen(n, BAM_CSOFT_CLIP);
      } else if (del_len < n && del_len > 0){
	ncigar[j] = bam_cigar_gen(del_len, BAM_CSOFT_CLIP);
      } else if (del_len == 0) {	// Adding insertions before start position of read
	ncigar[j] = bam_cigar_gen(n, BAM_CSOFT_CLIP);
	j++;
	i++;
	continue;
      }
      j++;
      ref_add = std::min(del_len, n);
      temp = n;
      n = std::max(n - del_len, 0);
      del_len = std::max(del_len - temp, 0);
      if(n > 0){
	ncigar[j] = bam_cigar_gen(n, cig);
	j++;
      }
      if(del_len ==0 && (bam_cigar_type(ncigar[j-1]) & 1) && (bam_cigar_type(ncigar[j-1]) & 2)){ // After soft clipping of query complete, keep incrementing start_pos until first base that consumes both query and ref
      	pos_start = true;
      }
    }
    if((bam_cigar_type(cig) & 2)) { // Consumes reference but not query
      start_pos += ref_add;
    }
    i++;
  }
  if(reverse){
    reverse_cigar(ncigar, j);
  }
  return {
    ncigar,
    true,
    j,
    start_pos
  };
}

void replace_cigar(bam1_t *b, uint32_t n, uint32_t *cigar){
  if (n != b->core.n_cigar) {
    int o = b->core.l_qname + b->core.n_cigar * 4;
    if (b->l_data + (n - b->core.n_cigar) * 4 > b->m_data) {
      b->m_data = b->l_data + (n - b->core.n_cigar) * 4;
      kroundup32(b->m_data);
      b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->l_data - o);
    memcpy(b->data + b->core.l_qname, cigar, n * 4);
    b->l_data += (n - b->core.n_cigar) * 4;
    b->core.n_cigar = n;
  } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}

void print_primers(std::vector<primer> primers){
  for(std::vector<primer>::iterator it = primers.begin(); it != primers.end(); ++it) {
    std::cout << "Get Start " << it->get_start() << "\n";
    std::cout << "Get End " << it->get_end() << "\n";
    std::cout << "Index " << it->get_indice() << "\n";
  }
}

int binarySearch(std::vector<primer> primers, uint32_t item, int low, int high){
  while (low <= high) {
    int mid = low + (high-low)/2;
    if (item == primers[mid].get_start()){
        return mid+1;
    } else if (item > primers[mid].get_start()){
        low = mid + 1;
    } else {
        high = mid -1;
    }
  }
  return low;
}

std::vector<primer> insertionSort(std::vector<primer> primers, uint32_t n){
  //std::cout << "in insertion sort length: " << n << "\n";
  uint32_t i = 0;
  int loc = 0;
  int j = 0;
 
  //iterate over vector of primers
  for (i=1; i < n; ++i){
    j = i-1;
    primer selected = primers[i];
    
    loc = binarySearch(primers, primers[i].get_start(), 0, j);
    while(j >= loc){
      primers[j+1] = primers[j];
      j--;
    }
    primers[j+1] = selected;
  }
  return primers;
}

// For paired reads
void get_overlapping_primers(bam1_t* r, std::vector<primer> primers, std::vector<primer> &overlapped_primers){
  overlapped_primers.clear();
  uint32_t start_pos = -1;
  char strand = '+';
  if(bam_is_rev(r)){
    start_pos = bam_endpos(r)-1;
    strand = '-';
  } else {
    start_pos = r->core.pos;
  }
   
  //print_primers(primers);
  //sort it first
  std::vector<primer> test = insertionSort(primers, primers.size());  
  //std::cout << test.size() << "\n";
  //then we iterate and push what fits
  for(std::vector<primer>::iterator it = test.begin(); it != test.end(); ++it) {
    //if we've passed the end, we're going to find no more matches
    if(start_pos < it->get_start()){
        //std::cout << "break start_pos: " << start_pos << " start:  " << it->get_end() << " end: " << it->get_start() << "\n";
        break;
    }
    if(start_pos >= it->get_start() && start_pos <= it->get_end() && (strand == it->get_strand() || it->get_strand() == 0))
      //std::cout << it->get_start() << " " << start_pos << "\n";
      overlapped_primers.push_back(*it);
  }
  //print_primers(overlapped_primers);
  //std::cout << "break\n";
}

// For unpaired reads
void get_overlapping_primers(bam1_t* r, std::vector<primer> primers, std::vector<primer> &overlapped_primers, bool unpaired_rev){
  overlapped_primers.clear();
  uint32_t start_pos = -1;
  char strand = '+';
  if(unpaired_rev){
    start_pos = bam_endpos(r) - 1;
    strand = '-';
  } else {
    start_pos = r->core.pos;
  }
  for(std::vector<primer>::iterator it = primers.begin(); it != primers.end(); ++it) {
    if(start_pos >= it->get_start() && start_pos <= it->get_end() && (strand == it->get_strand() ||it->get_strand() == 0))
      overlapped_primers.push_back(*it);
  }
}

void condense_cigar(cigar_ *t){
  uint32_t i = 0, len = 0, cig, next_cig;
  while(i< t->nlength -1){
    cig = bam_cigar_op(t->cigar[i]);
    next_cig = bam_cigar_op(t->cigar[i+1]);
    if(cig == next_cig){
      len = bam_cigar_oplen(t->cigar[i])+bam_cigar_oplen(t->cigar[i+1]);
      t->cigar[i] = bam_cigar_gen(len, bam_cigar_op(t->cigar[i]));
      for(uint32_t j = i+1; j < t->nlength - 1; j++){
	t->cigar[j] = t->cigar[j+1];
      }
      t->nlength--;
    } else {
      i++;
    }
  }
}

void add_pg_line_to_header(bam_hdr_t** hdr, char *cmd){
  size_t len = strlen((*hdr)->text) + strlen(cmd)+1;
  char * new_text = (char *)malloc(len);
  memcpy(new_text, (*hdr)->text, strlen((*hdr)->text));
  new_text[strlen((*hdr)->text)] = '\0';
  strcat(new_text, cmd);
  free((*hdr)->text);
  (*hdr)->text = new_text;
  new_text = NULL;
  (*hdr)->l_text = len-1;
}

// get the length of the longest primer
int get_bigger_primer(std::vector<primer> primers){
  int max_primer_len = 0;
  for (auto & p : primers) {
    if(max_primer_len < p.get_length()){
      max_primer_len = p.get_length();
    }
  }
  return max_primer_len;
}

// check if read is enveloped by any of the amplicons
bool amplicon_filter(IntervalTree amplicons, bam1_t* r){
  Interval fragment_coords = Interval(0, 1);
  if(r->core.isize > 0){
    fragment_coords.low = r->core.pos;
    fragment_coords.high = r->core.pos + r->core.isize;
  } else {
    fragment_coords.low = bam_endpos(r) + r->core.isize;
    fragment_coords.high = bam_endpos(r);
  }
  // debugging
  bool amplicon_flag = amplicons.envelopSearch(fragment_coords);
  return amplicon_flag;
}

int trim_bam_qual_primer(std::string bam, std::string bed, std::string bam_out, std::string region_, uint8_t min_qual, uint8_t sliding_window, std::string cmd, bool write_no_primer_reads, bool keep_for_reanalysis, int min_length = 30, std::string pair_info = "", int32_t primer_offset = 0) {
  int retval = 0;
  std::vector<primer> primers;
  int max_primer_len = 0;
  if(!bed.empty()){
    primers = populate_from_file(bed, primer_offset);
    if(primers.size() == 0){
      std::cout << "Exiting." << std::endl;
      return -1;
    }
  }
  max_primer_len = get_bigger_primer(primers);
  // get coordinates of each amplicon
  IntervalTree amplicons;
  if(!pair_info.empty()){
    amplicons = populate_amplicons(pair_info, primers);
    std::cout << "Amplicons detected: " << std::endl;
    amplicons.inOrder();
  }
  
  if(bam.empty()){
    std::cout << "Bam file is empty." << std::endl;
    return -1;
  }
  bam_out += ".bam";
  samFile *in = hts_open(bam.c_str(), "r");
  BGZF *out = bgzf_open(bam_out.c_str(), "w");
  if(in == NULL) {
    std::cout << ("Unable to open BAM file.") << std::endl;
    return -1;
  }
  //Load the index
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  if(idx == NULL) {
    std::cout << "Building BAM index" << std::endl;
    if(sam_index_build2(bam.c_str(), 0, 0)< 0){
      std::cout << ("Unable to open or build BAM index.") << std::endl;
      return -1;
    } else {
      idx = sam_index_load(in, bam.c_str());
    }
  }
  //Get the header
  bam_hdr_t *header = sam_hdr_read(in);
  if(header == NULL) {
    sam_close(in);
    std::cout << "Unable to open BAM header." << std::endl;
  }
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  if(bam_hdr_write(out, header) < 0){
    std::cout << "Unable to write BAM header to path." << std::endl;
    sam_close(in);
    return -1;
  }
  // Get relevant region
  int region_id = -1;
  uint64_t unmapped, mapped, log_skip;
  std::cout << std::endl << "Number of references in file: " << header->n_targets << std::endl;
  for (int i = 0; i < header->n_targets; ++i){
    std::cout << header->target_name[i] << std::endl;
    if(region_.compare(std::string(header->target_name[i])) == 0){
      region_id = i;
    }
    if(i==0){			// Reading only first reference
      region_.assign(header->target_name[i]);
      region_id = i;
    }
  }
  std::cout << "Using Region: " << region_ << std::endl << std::endl;
  // Get index stats
  hts_idx_get_stat(idx, region_id, &mapped, &unmapped);
  std::cout << "Found " << mapped << " mapped reads" << std::endl;
  std::cout << "Found " << unmapped << " unmapped reads" << std::endl;
  std::string hdr_text(header->text);
  if (hdr_text.find(std::string("SO:coordinate")) != std::string::npos) {
    std::cout << "Sorted By Coordinate" << std::endl; // Sort by coordinate
  } else if(hdr_text.find(std::string("SO:queryname")) != std::string::npos) {
    std::cout << "Sorted By Query Name" << std::endl; // Sort by name
  } else {
    std::cout << "Not sorted" << std::endl;
  }
  std::cout << "-------" << std::endl;
  log_skip = (mapped + unmapped > 10) ? (mapped + unmapped)/10 : 2;
  //Initialize iterator
  hts_itr_t *iter = NULL;
  //Move the iterator to the region we are interested in
  iter  = sam_itr_querys(idx, header, region_.c_str());
  if(header == NULL || iter == NULL) {
    sam_close(in);
    std::cout << "Unable to iterate to region within BAM/SAM." << std::endl;
    return -1;
  }
  //Initiate the alignment record
  bam1_t *aln = bam_init1();
  int ctr = 0;
  cigar_ t;
  init_cigar(&t);
  uint32_t primer_trim_count = 0, no_primer_counter = 0, low_quality = 0;
  bool unmapped_flag = false;
  bool amplicon_flag = false;
  bool isize_flag = true;
  uint32_t failed_frag_size = 0;
  uint32_t unmapped_counter = 0;
  uint32_t amplicon_flag_ctr = 0;
  primer cand_primer;
  std::vector<primer> overlapping_primers;
  std::vector<primer>::iterator cit;
  bool primer_trimmed = false;
  //Iterate through reads
  while(sam_itr_next(in, iter, aln) >= 0) {
    unmapped_flag = false;
    primer_trimmed = false;
    get_overlapping_primers(aln, primers, overlapping_primers);
    if((aln->core.flag&BAM_FUNMAP) == 0){ // If mapped
      // if primer pair info provided, check if read correctly overlaps with atleast one amplicon
      if(!pair_info.empty()){
        amplicon_flag = amplicon_filter(amplicons, aln);
        if(!amplicon_flag){
	  if (keep_for_reanalysis) {   // -k (keep) option
	    aln->core.flag |= BAM_FQCFAIL;
	    if (bam_write1(out, aln) < 0) { retval = -1; goto error; }
          }
          amplicon_flag_ctr++;
          continue;
	}
      }
      isize_flag = (abs(aln->core.isize) - max_primer_len) > abs(aln->core.l_qseq);
      // if reverse strand
      if((aln->core.flag&BAM_FPAIRED) != 0 && isize_flag){ // If paired
	get_overlapping_primers(aln, primers, overlapping_primers);
	if(overlapping_primers.size() > 0){ // If read starts before overlapping regions (?)
	  primer_trimmed = true;
	  if(bam_is_rev(aln)){	// Reverse read
	    cand_primer = get_min_start(overlapping_primers); // fetch reverse primer (?)
	    t = primer_trim(aln, isize_flag, cand_primer.get_start() - 1, false);
	  } else {		// Forward read
	    cand_primer = get_max_end(overlapping_primers); // fetch forward primer (?)
	    t = primer_trim(aln, isize_flag, cand_primer.get_end() + 1, false);
	    aln->core.pos += t.start_pos;
	  }
	  replace_cigar(aln, t.nlength, t.cigar);
	  free(t.cigar);
	  // Add count to primer
	  cit = std::find(primers.begin(), primers.end(), cand_primer);
	  if(cit != primers.end())
	    cit->add_read_count(1);
	}
	t = quality_trim(aln, min_qual, sliding_window);	// Quality Trimming
	if(bam_is_rev(aln))  // if reverse strand
	  aln->core.pos = t.start_pos;
	condense_cigar(&t);
	// aln->core.pos += t.start_pos;
	replace_cigar(aln, t.nlength, t.cigar);
      } else {			// Unpaired reads: Might be stitched reads
        if(abs(aln->core.isize) <= abs(aln->core.l_qseq)){
          failed_frag_size++;
        }
	// Forward primer
	get_overlapping_primers(aln, primers, overlapping_primers, false);
	if(overlapping_primers.size() > 0){
	  primer_trimmed = true;
	  cand_primer = get_max_end(overlapping_primers);
	  t = primer_trim(aln, isize_flag, cand_primer.get_end() + 1, false);
	  // Update read's left-most coordinate
	  aln->core.pos += t.start_pos;
	  replace_cigar(aln, t.nlength, t.cigar);
	  // Add count to primer
	  cit = std::find(primers.begin(), primers.end(), cand_primer);
	  if(cit != primers.end())
	    cit->add_read_count(1);
	}
	// Reverse primer
	get_overlapping_primers(aln, primers, overlapping_primers, true);
	if(overlapping_primers.size() > 0){
	  primer_trimmed = true;
	  cand_primer = get_min_start(overlapping_primers);
	  t = primer_trim(aln, isize_flag, cand_primer.get_start() - 1, true);
	  replace_cigar(aln, t.nlength, t.cigar);
	  // Add count to primer
	  cit = std::find(primers.begin(), primers.end(), cand_primer);
	  if(cit != primers.end())
	    cit->add_read_count(1);
	}
	t = quality_trim(aln, min_qual, sliding_window);	// Quality Trimming
	if(bam_is_rev(aln))  // if reverse strand
	  aln->core.pos = t.start_pos;
	condense_cigar(&t);
	replace_cigar(aln, t.nlength, t.cigar);
      }
      if(primer_trimmed){
	primer_trim_count++;
      }
    } else {
      unmapped_flag = true;
      unmapped_counter++;
      continue;
    }
    if(bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln)) >= min_length){
      if(primer_trimmed){	// Write to BAM only if primer found.
	int16_t cand_ind = cand_primer.get_indice();
	bam_aux_append(aln, "XA", 's', sizeof(cand_ind), (uint8_t*) &cand_ind);
	if(bam_write1(out, aln) < 0) { retval = -1; goto error; }
      } else {  // no primer found
        if (keep_for_reanalysis) {   // -k (keep) option
          if((primers.size() == 0 || !write_no_primer_reads) && !unmapped_flag){ // -k only option
            aln->core.flag |= BAM_FQCFAIL;
          }
	  if (bam_write1(out, aln) < 0) { retval = -1; goto error; }
	} else {        // no -k option
          if((primers.size() == 0 || write_no_primer_reads) && !unmapped_flag){ // -e only option
	    if (bam_write1(out, aln) < 0) { retval = -1; goto error; }
          }
        }
	no_primer_counter++;
      }
    } else {
      low_quality++;
      if (keep_for_reanalysis) {
        aln->core.flag |= BAM_FQCFAIL;
	if (bam_write1(out, aln) < 0) { retval = -1; goto error; }
      }
    }
    ctr++;
    if(ctr % log_skip == 0){
      std::cout << "Processed " << (ctr/log_skip) * 10 << "% reads ... " << std::endl;
    }
  }
  std::cout << std::endl << "-------" << std::endl;
  std::cout << "Results: " << std::endl;
  std::cout << "Primer Name" << "\t" << "Read Count" << std::endl;
  for(cit = primers.begin(); cit != primers.end(); ++cit) {
    std::cout << cit->get_name() << "\t" << cit->get_read_count() << std::endl;
  }
  std::cout << std::endl << "Trimmed primers from " << round_int(primer_trim_count, mapped) << "% (" << primer_trim_count <<  ") of reads." << std::endl;
  std::cout << round_int( low_quality, mapped) << "% (" << low_quality << ") of reads were quality trimmed below the minimum length of " << min_length << " bp and were ";
  if (keep_for_reanalysis) {
    std::cout << "marked as failed" << std::endl;
  } else {
    std::cout << "not written to file." << std::endl;
  }
  if(write_no_primer_reads){
    std::cout << round_int(no_primer_counter, mapped) << "% ("  << no_primer_counter << ")"
              << " of reads started outside of primer regions. Since the "
              << (keep_for_reanalysis ? "-ek flags were " : "-e flag was ")
              << "given, these reads were written to file";
    std::cout << "." << std::endl;
  } else if (primers.size() == 0) {
    std::cout << round_int(no_primer_counter, mapped) << "% ("  << no_primer_counter << ") of reads started outside of primer regions. Since there were no primers found in BED file, these reads were written to file." << std::endl;
  } else {
    std::cout << round_int(no_primer_counter, mapped) << "% ("  << no_primer_counter
              << ") of reads that started outside of primer regions were ";
    if (keep_for_reanalysis) {
      std::cout << "written to file and marked as failed";
    } else {
      std::cout << "not written to file";
    }
    std::cout << std::endl;
  }
  if(unmapped_counter > 0){
    std::cout << unmapped_counter << " unmapped reads were not written to file." << std::endl;
  }
  if(amplicon_flag_ctr > 0){
    std::cout << round_int(amplicon_flag_ctr, mapped) 
              << "% (" << amplicon_flag_ctr 
              << ") reads were ignored because they did not fall within an amplicon" 
              << std::endl;
  }
  if(failed_frag_size > 0){
    std::cout << round_int(failed_frag_size, mapped)
              << "% (" << failed_frag_size
              << ") of reads had their insert size smaller than their read length"
              << std::endl;
  }

 error:
  if (retval) std::cout << "Not able to write to BAM" << std::endl;
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  bgzf_close(out);
  return retval;
}
