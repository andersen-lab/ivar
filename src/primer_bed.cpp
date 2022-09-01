#include "primer_bed.h"

std::string primer::get_name(){
  return name;
}

std::string primer::get_region(){
  return region;
}

int primer::get_score(){
  return score;
}

uint32_t primer::get_start() const {
  return start;
}

uint32_t primer::get_end() const{
  return end;
}

char primer::get_strand(){
  return strand;
}

int primer::get_length(){
  return end - start + 1;
}

int16_t primer::get_pair_indice(){
  return pair_indice;
}

int16_t primer::get_indice() const{
  return indice;
}

uint32_t primer::get_read_count() const{
  return read_count;
}

void primer::set_start(uint32_t  s){
  start = s;
}

void primer::set_end(uint32_t e){
  end = e;
}

void primer::set_strand(char s){
  strand = s;
}

void primer::set_region(std::string r){
  region = r;
}

void primer::set_name(std::string n){
  name = n;
}

void primer::set_score(int s){
  score = s;
}

void primer::set_pair_indice(int16_t i){
  pair_indice = i;
}

void primer::set_indice(int16_t i){
  indice = i;
}

void primer::set_read_count(uint32_t rc){
  read_count = rc;
}

void primer::add_read_count(uint32_t rc){
  read_count += rc;
}

void print_bed_format(){
  std::cout << "iVar uses the standard 6 column BED format as defined here - https://genome.ucsc.edu/FAQ/FAQformat.html#format1." << std::endl;
  std::cout << "It requires the following columns delimited by a tab: chrom, chromStart, chromEnd, name, score, strand" << std::endl;
}

std::vector<primer> populate_from_file(std::string path, int32_t offset = 0){
  std::ifstream  data(path.c_str());
  std::string line;
  std::vector<primer> primers;
  int16_t indice = 0;
  while(std::getline(data,line)){ // Remove extra lineStream
    std::stringstream lineStream(line);
    std::string cell;
    int ctr = 0;
    primer p;
    p.set_strand(0);		// Set strand to NULL
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	p.set_region(cell);
	break;
      case 1:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_start(std::stoul(cell) - offset);
	} else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
	break;
      case 2:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_end(std::stoul(cell) - 1 + offset); // Bed format - End is not 0 based
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
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_score(stoi(cell));
	} else {
	  print_bed_format();  // score is missing, send warning but continue populating
    std::cout << "\nWARNING: The BED file provided did not have the expected score column, but iVar will continue trimming\n" << std::endl;
    p.set_score(-1);
	}
	break;
      case 5:
	if(cell[0] == '+' || cell[0] == '-')
	  p.set_strand(cell[0]);
	else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
      }
      ctr++;
    }
    if(indice == 0 && ctr < 6)
      std::cout << "Strand not found in primer BED file so strand will not be considered for trimming" << std::endl;
    p.set_indice(indice);
    p.set_pair_indice(-1);
    p.set_read_count(0);
    primers.push_back(p);
    indice++;
  }
  std::cerr << "Found " << primers.size() << " primers in BED file" << std::endl;
  return primers;
}

std::vector<primer> populate_from_file(std::string path){
  std::ifstream  data(path.c_str());
  std::string line;
  std::vector<primer> primers;
  int16_t indice = 0;
  while(std::getline(data,line)){ // Remove extra lineStream
    std::stringstream lineStream(line);
    std::string cell;
    int ctr = 0;
    primer p;
    p.set_strand(0);		// Set strand to NULL
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	p.set_region(cell);
	break;
      case 1:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_start(std::stoul(cell));
	} else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
	break;
      case 2:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_end(std::stoul(cell) - 1); // Bed format - End is not 0 based
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
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_score(stoi(cell));
	} else {
	  print_bed_format();  // score is missing, send warning but continue populating
    std::cout << "\nWARNING: The BED file provided did not have the expected score column, but iVar will continue trimming\n" << std::endl;
    p.set_score(-1);
	}
	break;
      case 5:
	if(cell[0] == '+' || cell[0] == '-')
	  p.set_strand(cell[0]);
	else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
      }
      ctr++;
    }
    if(indice == 0 && ctr < 6)
      std::cout << "Strand not found in primer BED file so strand will not be considered for trimming" << std::endl;
    p.set_indice(indice);
    p.set_pair_indice(-1);
    p.set_read_count(0);
    primers.push_back(p);
    indice++;
  }
  std::cout << "Found " << primers.size() << " primers in BED file" << std::endl;
  return primers;
}

std::vector<primer> get_primers(std::vector<primer> p, unsigned int pos){
  std::vector<primer> primers_with_mismatches;
  for(std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    if(it->get_start() <= pos && it->get_end() >= pos){
      primers_with_mismatches.push_back(*it);
    }
  }
  return primers_with_mismatches;
}

// Assumes unique primer names in BED file
int get_primer_indice(std::vector<primer> p, std::string name){
  for(std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    if(it->get_name().compare(name) == 0){
      return it - p.begin();
    }
  }
  return -1;
}

int populate_pair_indices(std::vector<primer> &primers, std::string path){
  std::ifstream fin(path.c_str());
  std::string line, cell, p1,p2;
  std::stringstream line_stream;
  std::vector<primer>::iterator it;
  int32_t indice;
  while (std::getline(fin, line)){
    line_stream << line;
    std::getline(line_stream, cell, '\t');
    p1 = cell;
    line_stream.clear();
    std::getline(line_stream, cell, '\t');
    p2 = cell;
    line_stream.clear();
    if(!p1.empty() && !p2.empty()){
      for(it = primers.begin(); it != primers.end(); ++it) {
	if (it->get_name() == p1) {
	  indice = get_primer_indice(primers, p2);
	  if (indice != -1)
	    it->set_pair_indice(indice);
	  else
	    std::cout << "Primer pair for " << p1 << " not found in BED file." <<std::endl;
	} else if (it->get_name() == p2){
	  indice = get_primer_indice(primers, p1);
	  if(indice != -1)
	    it->set_pair_indice(indice);
	  else
	    std::cout << "Primer pair for " << p2 << " not found in BED file." << std::endl;
	}
      }
    }
  }
  return 0;
}

primer get_min_start(std::vector<primer> primers){
  std::vector<primer>::iterator it;
  auto minmax_start = std::minmax_element(primers.begin(), primers.end(), [] (primer lhs, primer rhs) {return lhs.get_start() < rhs.get_start();});
  return *(minmax_start.first);
}

primer get_max_end(std::vector<primer> primers){
  auto minmax_start = std::minmax_element(primers.begin(), primers.end(), [] (primer lhs, primer rhs) {return lhs.get_end() < rhs.get_end();});
  return *(minmax_start.second);
}