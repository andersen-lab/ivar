#include "primer_bed.h"
#include "allele_functions.h"
#include <chrono>
std::string primer::get_name() { return name; }

std::string primer::get_region() { return region; }

std::vector<std::vector<uint32_t>> cigarotype::get_cigarotypes() { return cigarotypes; }

std::vector<std::vector<uint8_t>> cigarotype::get_aux_tags() { return aux_tags; }

std::vector<std::vector<uint32_t>> cigarotype::get_qualities() { return qualities; }

std::vector<std::vector<uint8_t>> cigarotype::get_sequences() { return sequences; }

std::vector<std::string> cigarotype::get_qnames() { return qnames; }

std::vector<uint32_t> cigarotype::get_count_cigarotypes() { return count_cigarotypes; }

std::vector<uint32_t> cigarotype::get_start_positions() { return ncigarotypes; }

std::vector<uint32_t> cigarotype::get_nlengths() { return nlengths; }

int primer::get_score() { return score; }

std::vector<position> primer::get_positions() { return positions; }

void primer::set_positions(position pos) { positions.push_back(pos); }

uint32_t primer::get_start() const { return start; }

uint32_t primer::get_end() const { return end; }

char primer::get_strand() { return strand; }

int primer::get_length() { return end - start + 1; }

int16_t primer::get_pair_indice() { return pair_indice; }

int16_t primer::get_indice() const { return indice; }

uint32_t primer::get_read_count() const { return read_count; }

void primer::set_start(uint32_t s) { start = s; }

void primer::set_end(uint32_t e) { end = e; }

void primer::set_strand(char s) { strand = s; }

void primer::set_region(std::string r) { region = r; }

void primer::set_name(std::string n) { name = n; }

void primer::set_score(int s) { score = s; }

void primer::set_pair_indice(int16_t i) { pair_indice = i; }

void primer::set_indice(int16_t i) { indice = i; }

void primer::set_read_count(uint32_t rc) { read_count = rc; }

void primer::add_read_count(uint32_t rc) { read_count += rc; }

void print_bed_format() {
  std::cerr << "iVar uses the standard 6 column BED format as defined here - "
               "https://genome.ucsc.edu/FAQ/FAQformat.html#format1."
            << std::endl;
  std::cerr << "It requires the following columns delimited by a tab: chrom, "
               "chromStart, chromEnd, name, score, strand"
            << std::endl;
}

void primer::populate_positions(){
  for(uint32_t i=0; i < this->hardcoded_length; i++){
    position tmp;
    tmp.alleles = populate_basic_alleles();
    if(this->get_strand() == '+'){
      tmp.pos = i + this->get_start();
    } else {
      if(this->get_end() > i){
        break;
      }
      tmp.pos = this->get_end() - i;
    }
    positions.push_back(tmp);
  }
}

std::vector<primer> populate_from_file(std::string path, int32_t offset = 0) {
  std::ifstream data(path.c_str());
  std::string line;
  std::vector<primer> primers;
  int16_t indice = 0;
  while (std::getline(data, line)) {  // Remove extra lineStream
    std::stringstream lineStream(line);
    std::string cell;
    int ctr = 0;
    primer p;
    p.set_strand(0);  // Set strand to NULL
    while (std::getline(lineStream, cell, '\t')) {
      switch (ctr) {
        case 0:
          p.set_region(cell);
          break;
        case 1:
          if (std::all_of(cell.begin(), cell.end(), ::isdigit)) {
            p.set_start(std::stoul(cell) - offset);
          } else {
            print_bed_format();
            primers.clear();
            return primers;
          }
          break;
        case 2:
          if (std::all_of(cell.begin(), cell.end(), ::isdigit)) {
            p.set_end(std::stoul(cell) - 1 +
                      offset);  // Bed format - End is not 0 based
          } else {
            print_bed_format();
            primers.clear();
            return primers;
          }
          break;
        case 3:
          p.set_name(cell);
          break;
        case 4:
          if (std::all_of(cell.begin(), cell.end(), ::isdigit)) {
            p.set_score(stoi(cell));
          } else {
            print_bed_format();  // score is missing, send warning but continue
                                 // populating
            std::cerr
                << "\nWARNING: The BED file provided did not have the expected "
                   "score column, but iVar will continue trimming\n"
                << std::endl;
            p.set_score(-1);
          }
          break;
        case 5:
          if (cell[0] == '+' || cell[0] == '-')
            p.set_strand(cell[0]);
          else {
            print_bed_format();
            primers.clear();
            return primers;
          }
      }
      ctr++;
    }
    if (indice == 0 && ctr < 6)
      std::cerr << "Strand not found in primer BED file so strand will not be "
                   "considered for trimming"
                << std::endl;
    p.set_indice(indice);
    p.set_pair_indice(-1);
    p.set_read_count(0);
    primers.push_back(p);
    indice++;
  }
  std::cerr << "Found " << primers.size() << " primers in BED file"
            << std::endl;
  return primers;
}

std::vector<primer> populate_from_file(std::string path) {
  std::ifstream data(path.c_str());
  std::string line;
  std::vector<primer> primers;
  int16_t indice = 0;
  while (std::getline(data, line)) {  // Remove extra lineStream
    std::stringstream lineStream(line);
    std::string cell;
    int ctr = 0;
    primer p;
    p.set_strand(0);  // Set strand to NULL
    while (std::getline(lineStream, cell, '\t')) {
      switch (ctr) {
        case 0:
          p.set_region(cell);
          break;
        case 1:
          if (std::all_of(cell.begin(), cell.end(), ::isdigit)) {
            p.set_start(std::stoul(cell));
          } else {
            print_bed_format();
            primers.clear();
            return primers;
          }
          break;
        case 2:
          if (std::all_of(cell.begin(), cell.end(), ::isdigit)) {
            p.set_end(std::stoul(cell) - 1);  // Bed format - End is not 0 based
          } else {
            print_bed_format();
            primers.clear();
            return primers;
          }
          break;
        case 3:
          p.set_name(cell);
          break;
        case 4:
          if (std::all_of(cell.begin(), cell.end(), ::isdigit)) {
            p.set_score(stoi(cell));
          } else {
            print_bed_format();  // score is missing, send warning but continue
                                 // populating
            std::cerr
                << "\nWARNING: The BED file provided did not have the expected "
                   "score column, but iVar will continue trimming\n"
                << std::endl;
            p.set_score(-1);
          }
          break;
        case 5:
          if (cell[0] == '+' || cell[0] == '-')
            p.set_strand(cell[0]);
          else {
            print_bed_format();
            primers.clear();
            return primers;
          }
      }
      ctr++;
    }
    if (indice == 0 && ctr < 6)
      std::cerr << "Strand not found in primer BED file so strand will not be "
                   "considered for trimming"
                << std::endl;
    p.set_indice(indice);
    p.set_pair_indice(-1);
    p.set_read_count(0);
    primers.push_back(p);
    indice++;
  }
  std::cerr << "Found " << primers.size() << " primers in BED file"
            << std::endl;
  return primers;
}

std::vector<primer> get_primers(std::vector<primer> p, unsigned int pos) {
  std::vector<primer> primers_with_mismatches;
  for (std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    if (it->get_start() <= pos && it->get_end() >= pos) {
      primers_with_mismatches.push_back(*it);
    }
  }
  return primers_with_mismatches;
}

// function to trim trailing spaces
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ") {
  str.erase(str.find_last_not_of(chars) + 1);
  return str;
}

// Assumes unique primer names in BED file
// returns the index of a primer by name
int get_primer_indice(std::vector<primer> p, std::string name) {
  // iterate through the primers from the bed file
  for (std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    // check if the two strings are the same
    if (it->get_name().compare(rtrim(name)) == 0) {
      return it - p.begin();
    }
  }
  return -1;
}

// function using the tab seperated primer pair file
int populate_pair_indices(std::vector<primer>& primers, std::string path) {
  /*
   * @param primers: the primer vector to add pair info to
   * @param path: the path to the primer pair file
   */

  // load primer pair file
  std::ifstream fin(path.c_str());
  std::string line, cell, p1, p2;
  std::stringstream line_stream;
  std::vector<primer>::iterator it;
  int32_t indice;
  // iterate the primer pair file line by line
  while (std::getline(fin, line)) {
    line_stream << line;
    std::getline(line_stream, cell, '\t');
    p1 = cell;
    line_stream.clear();
    std::getline(line_stream, cell, '\t');
    p2 = cell;
    line_stream.clear();

    p1 = rtrim(p1);
    p2 = rtrim(p2);

    if (!p1.empty() && !p2.empty()) {
      for (it = primers.begin(); it != primers.end(); ++it) {
        // search for primer name in pair file
        if (it->get_name() == p1) {
          // make sure it's pair exists
          indice = get_primer_indice(primers, p2);
          if (indice != -1) {
            it->set_pair_indice(indice);
          } else {
            std::cerr << "Primer pair for " << p1 << " not found in BED file."
                      << std::endl;
          }
        } else if (it->get_name() == p2) {
          indice = get_primer_indice(primers, p1);
          if (indice != -1)
            it->set_pair_indice(indice);
          else
            std::cerr << "Primer pair for " << p2 << " not found in BED file."
                      << std::endl;
        }
      }
    } else {
      std::cerr << "Primer pair is empty." << std::endl;
    }
  }
  return 0;
}

void primer::transform_mutations() {
  /*
   * Take all recorded cigarotypes and transform them into positions relative to primer
   */
  std::vector<std::vector<uint32_t>> cigarotypes = this->get_cigarotypes(); 
  std::vector<uint32_t> start_positions = this->get_start_positions();
  std::vector<std::vector<uint8_t>> sequences = this->get_sequences();
  std::vector<std::vector<uint8_t>> aux_tags = this->get_aux_tags();
  std::vector<uint32_t> counts = this->get_count_cigarotypes();
  std::vector<std::string> qnames = this->get_qnames();
  std::vector<std::vector<uint32_t>> qualities = this->get_qualities(); 
  //this tracks all mutations at all positions
  std::string test = "";
  //here let's turn the cigar string into a vector of alleles specific to this primer
  //iterate all unique sequences
  for(uint32_t i=0; i < cigarotypes.size(); i++){
    //get the cigar string and start pos
    std::vector<uint32_t> quality = qualities[i];
    std::vector<uint32_t> cigarotype = cigarotypes[i]; //carries the insertions
    std::vector<uint8_t> aux_tag = aux_tags[i]; //carries deletions
    std::vector<uint8_t> sequence = sequences[i]; //carries NT values
    std::vector<uint32_t> sc_positions; //sc positions            
    uint32_t ccount = counts[i];                                                     
    uint32_t start_pos = start_positions[i]; // pos after soft-clipped region
    std::string qname = qnames[i];                                             
    uint32_t consumed_query = 0;
    uint32_t consumed_ref = 0;
    std::string nuc;
    std::vector<uint32_t> ignore_sequence; //positions in query that are insertions
    //we'll use the cigar string to handle insertions
    for(uint32_t j=0; j < cigarotype.size(); j++){
      uint32_t op = bam_cigar_op(cigarotype[j]);
      uint32_t oplen = bam_cigar_oplen(cigarotype[j]);
      if (qname == test){
        std::cerr << op << " " << oplen << std::endl;
      }
      //honestly this whole bit could be better - more general
      //insertions
      if(op == 1){
        std::ostringstream convert;
        uint32_t q = 0;
        for(uint32_t k=0; k < oplen; k++) {
          //convert data type to get the characters
          ignore_sequence.push_back(k+consumed_ref+start_pos);
          convert << sequence[k+consumed_query];
          q += quality[k+consumed_query];
        }
        nuc.clear();
        nuc = convert.str();
        int avg_q = (int)q/nuc.size();
        std::string nuc = "+" + convert.str();
        //check if this position exists
        int exists = check_position_exists(start_pos+consumed_ref, positions);
        if (exists != -1) {
          positions[exists].update_alleles(nuc, ccount, avg_q, false);  
        } else {
          //add position to vector
          position add_pos;
          add_pos.pos = start_pos+consumed_ref; //don't add a position
          add_pos.update_alleles(nuc, ccount, avg_q, false);            
          positions.push_back(add_pos);
        }
        consumed_query += oplen;
        continue;        
      }
      if (!(bam_cigar_type(op) & 1) && !(bam_cigar_type(op) & 2)){
        continue;
      }
      //if we don't consume both query and reference
      if (!(bam_cigar_type(op) & 1) || !(bam_cigar_type(op) & 2)){
        if (op != 2){
          for(uint32_t k=0; k < oplen; k++) {
            //convert data type to get the characters
            sc_positions.push_back(k+consumed_query);
          }
        }
      }
      //consumes query
      if (bam_cigar_type(op) & 1){
        consumed_query += oplen;
      } 
      //consumes ref
      if (bam_cigar_type(op) & 2){
        consumed_ref += oplen;
      } 
    }
    //we will use the aux tag to handle deletions and substituitons
    bool deletion = false;
    bool substitution = false;
    uint32_t current_pos = start_pos;
    std::vector<uint32_t> deletion_positions; //the aux tag does NOT recognize insertions
    std::vector<uint32_t> substitutions; //handle substitutions                                              
    std::string gather_digits = "";
    std::string deleted_char;   
    std::string sub_char;
    uint32_t last_char = 0;
    for(uint32_t j=1; j < aux_tag.size(); j++){
      char character = (char) aux_tag[j];
      if(test == qname){
        std::cerr << "aux " << aux_tag[j] << " digits " << gather_digits << " current pos " << current_pos << " deletion " << deletion << " is digit " << isdigit(character) <<  std::endl;
      }
      if (character == '^'){
        current_pos += std::stoi(gather_digits);
        gather_digits = "";
        deletion = true;
      } else if (isdigit(character) && deletion) {
        deleted_char = "";
        deletion = false;
        gather_digits += character;
      } else if (isalpha(character) && deletion) {
        int exists = check_position_exists(current_pos, positions);
        if (exists != -1) {
          positions[exists].update_alleles("-", ccount, 0, false);  
        } else {
          //add position to vector
          position add_pos;
          add_pos.pos = current_pos; //don't add a position
          add_pos.update_alleles("-", ccount, 0, false);            
          positions.push_back(add_pos);
        } 
        deletion_positions.push_back(current_pos);
        current_pos += 1;
        deleted_char += character;
        deletion = true;    
      } else if (isdigit(character) && !deletion) {
        if(substitution){
          for(uint32_t z=0; z < sub_char.size(); z++){
            substitutions.push_back(current_pos + z + 1);
          }
          substitution = false;
          sub_char.clear();
        }
        if(last_char > 0){
          current_pos += last_char;
          last_char = 0;
        }  
        gather_digits += character;
      } else if (isalpha(character) && !deletion) {
        last_char += 1;
        substitution = true;
        sub_char += character;
        if(gather_digits.size() > 0){
          current_pos += std::stoi(gather_digits);
        }
        gather_digits = "";
      } 
    }
    //now that we know where the insertions and deletions are, let's just iterate the query sequence and add it in, skipping problem positions
    current_pos = start_pos;
    bool prev_insertion = false;
    std::ostringstream convert;
    nuc.clear();
    std::vector<uint32_t> seen_insertions;
    std::vector<uint32_t>::iterator i_it;

    //j is relative to the sequence and current pos to the reference
    for(uint32_t j=0; j < sequence.size(); j++){
      //auto start = std::chrono::high_resolution_clock::now();
      std::vector<uint32_t>::iterator it = find(deletion_positions.begin(), deletion_positions.end(), current_pos); 
      if (it != deletion_positions.end()) {
        current_pos += 1;
        j -= 1;
        continue;
      }
      it = find(ignore_sequence.begin(), ignore_sequence.end(), current_pos);
      i_it = find(seen_insertions.begin(), seen_insertions.end(), current_pos);
      
      //handle insertions
      if (it != ignore_sequence.end() && i_it == seen_insertions.end()) {
        if(qname == test){
          std::cerr << j << " " << current_pos << " insertion" << std::endl;
        }
        std::ostringstream convert;
        convert << sequence[j];
        nuc += convert.str();
        seen_insertions.push_back(current_pos);
        current_pos += 1;
        prev_insertion = true;
        continue;
      } else if (it == ignore_sequence.end()){
        //insertion is over 
        if(prev_insertion) {
          current_pos -= nuc.size();
        }
        prev_insertion = false;
        nuc = "";
      }
      //auto end = std::chrono::high_resolution_clock::now();
      //std::chrono::duration<double> duration = end - start;
      //std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;
      it = find(sc_positions.begin(), sc_positions.end(), j);
      if (it != sc_positions.end()){
        continue;
      }
      //auto end2 = std::chrono::high_resolution_clock::now();
      //std::chrono::duration<double> duration2 = end2 - end;
      //std::cout << "Execution time 2: " << duration2.count() << " seconds." << std::endl;
      current_pos += 1;
      std::ostringstream convert;
      bool ref = false;
      convert << sequence[j];
      std::string nuc = convert.str(); 
      if (std::find(substitutions.begin(), substitutions.end(), current_pos) == substitutions.end()){
        ref = true;
      }
      /*if(current_pos == 28912 && nuc == "T"){
        std::cerr << qname << " " << ref << std::endl;
      }*/
      int exists = check_position_exists(current_pos, positions);
      if (exists != -1) {
        positions[exists].update_alleles(nuc, ccount, quality[j], ref);  
      } else {
        position add_pos;
        add_pos.pos = current_pos; //don't add a position
        add_pos.update_alleles(nuc, ccount, quality[j], ref);            
        positions.push_back(add_pos);
      }
      //auto end3 = std::chrono::high_resolution_clock::now();
      //std::chrono::duration<double> duration3 = end3 - end2;
      //std::cout << "Execution time 3: " << duration3.count() << " seconds." << std::endl; 
    }
  }
}

//cigar must be passed by pointer due to array decay
void cigarotype::add_cigarotype(uint32_t *cigar , uint32_t start_pos, uint32_t nlength, uint8_t *seq, uint8_t *aux, std::string qname, uint8_t *quality) {
  bool found = false; //have we seen this before
  std::vector<uint8_t> saved_aux; //placeholder for saved sequence
  std::vector<uint32_t> saved_cigarotype; //placeholder for saved cigarotypes
  std::vector<uint8_t> saved_seq;
  uint32_t sp; //placeholder for saved starting position
  std::vector<uint32_t> qual_reformat;
  std::vector<uint32_t> cigar_reformat;
  std::vector<uint8_t> seq_reformat;
  std::vector<uint8_t> aux_reformat;
  std::string test = "";
  uint32_t length=0;
  uint32_t ll=0;
  std::vector<uint32_t> useful;

  //first handle the array decay aspect
  for(uint32_t i=0; i < nlength; i++){
    cigar_reformat.push_back(cigar[i]);
    uint32_t op = bam_cigar_op(cigar[i]);
    //consumes query
    if (bam_cigar_type(op) & 1){
      length = bam_cigar_oplen(cigar[i]);
      for(uint32_t k=0; k < length; k++){
        useful.push_back(ll+k);
      }
      ll += length;
    } 
  }
  int i = 0;
  bool substitution = false;
  do{
    aux_reformat.push_back(aux[i]); 
    if(isalpha(aux[i]) && aux[i] != 'Z'){
      substitution = true;
    }
    i++;
  } while(aux[i] != '\0');

  for(uint32_t k=0; k < useful.size(); k++){
    qual_reformat.push_back(quality[useful[k]]);
    if(substitution){
      uint8_t nt = 0;
      nt = seq_nt16_str[bam_seqi(seq, useful[k])]; //this operation is costly....
      seq_reformat.push_back(nt);
    }
  }
  //cigarotypes have insertions
  //aux tags have deletions and location of substitution
  //nt seq had the actual value of substitution
  for(uint32_t i=0; i < cigarotypes.size(); i++){
    sp = ncigarotypes[i]; 
    saved_aux = aux_tags[i];
    saved_cigarotype = cigarotypes[i];
    saved_seq = sequences[i];
    if(cigar_reformat == saved_cigarotype && start_pos == sp && aux_reformat == saved_aux){
      if(substitution && seq_reformat != saved_seq){
        continue;
      }
      found = true;
      count_cigarotypes[i] += 1;
      for(uint32_t k = 0; k < qual_reformat.size(); k++){
        qualities[i][k] += (uint32_t) qual_reformat[k];
      }
      break;
    }
  }
  //haven't seen this cigar/start pos combo before
  if(!found){
    if(seq_reformat.size() == 0){ 
      uint8_t nt = 0;
      for(uint32_t k=0; k < useful.size(); k++){
        nt = seq_nt16_str[bam_seqi(seq, useful[k])];
        seq_reformat.push_back(nt);
      }
    }
    qualities.push_back(qual_reformat);
    cigarotypes.push_back(cigar_reformat);    
    ncigarotypes.push_back(start_pos);
    sequences.push_back(seq_reformat);
    aux_tags.push_back(aux_reformat);
    qnames.push_back(qname);

    uint32_t digit = 1;
    count_cigarotypes.push_back(digit);
    nlengths.push_back(nlength);
  }
}

primer get_min_start(std::vector<primer> primers) {
  std::vector<primer>::iterator it;
  auto minmax_start = std::minmax_element(
      primers.begin(), primers.end(),
      [](primer lhs, primer rhs) { return lhs.get_start() < rhs.get_start(); });
  return *(minmax_start.first);
}

primer get_max_end(std::vector<primer> primers) {
  auto minmax_start = std::minmax_element(
      primers.begin(), primers.end(),
      [](primer lhs, primer rhs) { return lhs.get_end() < rhs.get_end(); });
  return *(minmax_start.second);
}

void print_all_primer_info(std::vector<primer> primers) {
  std::vector<primer>::iterator it;
  for (it = primers.begin(); it != primers.end(); ++it) {
    std::cerr << "Primer name: " << it->get_name() << std::endl;
    std::cerr << "Primer start: " << it->get_start() << std::endl;
    std::cerr << "Primer end: " << it->get_end() << std::endl;
    std::cerr << "Indice: " << it->get_indice() << std::endl;
    std::cerr << "Pair indice: " << it->get_pair_indice() << std::endl;
  }
}

void print_primer_info(primer primer) {
  std::cerr << "Primer name: " << primer.get_name() << std::endl;
  std::cerr << "Primer start: " << primer.get_start() << std::endl;
  std::cerr << "Primer end: " << primer.get_end() << std::endl;
  std::cerr << "Indice: " << primer.get_indice() << std::endl;
  std::cerr << "Pair indice: " << primer.get_pair_indice() << std::endl;
}
