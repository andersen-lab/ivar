#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "htslib/sam.h"
#include "allele_functions.h"
#ifndef primer_bed
#define primer_bed

class cigarotype{
  protected:
    std::vector<std::vector<uint32_t>> cigarotypes; //unique cigar value 
    std::vector<uint32_t> nlengths; //this is just for printing, record length og cig ops TODO can remove
    std::vector<uint32_t> ncigarotypes; //starting pos of these cigars
    std::vector<uint32_t> count_cigarotypes; //count of amount found
    std::vector<std::vector<uint8_t>> aux_tags; //auxillary tags
    std::vector<std::vector<uint8_t>> sequences; //sequences                        
    std::vector<std::vector<uint32_t>> qualities; //qualities                                             
    std::vector<std::vector<std::string>> qnames; //TODO remove this at the end, qname
  public:
    std::vector<std::vector<uint32_t>> get_cigarotypes();
    std::vector<std::vector<uint8_t>> get_aux_tags();
    std::vector<std::vector<uint8_t>> get_sequences();   
    std::vector<std::vector<uint32_t>> get_qualities();
    std::vector<uint32_t> get_count_cigarotypes();
    std::vector<uint32_t> get_start_positions();
    std::vector<uint32_t> get_nlengths();
    std::vector<std::vector<std::string>> get_qnames(); 
    void add_cigarotype(uint32_t *cigar, uint32_t start_pos, uint32_t nlength, uint8_t *seq, uint8_t *aux, std::string qname, uint8_t *quality);
    int test = 12;

};

class primer : public cigarotype{
 private:
  std::string region;
  uint32_t start;  // 0 based
  uint32_t end;    // 0 based
  std::string name;
  int score = 0;
  char strand;
  int16_t pair_indice;
  int16_t indice;
  uint32_t read_count = 0;
  std::vector<position> positions; //allele counts for this primer

 public:
  uint32_t hardcoded_length=300;
  std::string get_name();
  std::string get_region();
  std::vector<position> get_positions();
  void set_positions(position pos);
  int get_score();
  uint32_t get_start() const;
  uint32_t get_end() const;
  char get_strand();
  int get_length();
  int16_t get_pair_indice();
  int16_t get_indice() const;
  uint32_t get_read_count() const;
  void populate_positions();
  void set_start(uint32_t s);
  void set_end(uint32_t e);
  void set_strand(char s);
  void set_region(std::string r);
  void set_name(std::string n);
  void set_score(int s);
  void set_pair_indice(int16_t i);
  void set_indice(int16_t i);
  void set_read_count(uint32_t rc);
  void add_read_count(uint32_t rc);
  void transform_mutations();
  bool operator==(const primer& p) const {
    return (indice == p.get_indice()) ? true : false;
  }
};

std::vector<primer> populate_from_file(std::string path, int32_t offset);
std::vector<primer> populate_from_file(std::string path);
std::vector<primer> get_primers(std::vector<primer> p, unsigned int pos);
int get_primer_indice(std::vector<primer> p, std::string name);
int populate_pair_indices(std::vector<primer>& primers, std::string path);
void print_primer_info(primer primers);
void print_all_primer_info(std::vector<primer> primers);
primer get_min_start(std::vector<primer> primers);
primer get_max_end(std::vector<primer> primers);
void add_cigarotype(uint32_t *cigar, uint32_t start_pos, uint32_t nlength, uint8_t *seq, uint8_t *aux, std::string qname, uint8_t *quality);
std::vector<std::vector<uint32_t>> get_cigarotypes();
std::vector<std::vector<uint8_t>> get_aux_tags();
std::vector<std::vector<uint32_t>> get_qualities();   
std::vector<std::vector<uint8_t>> get_sequences();   
std::vector<uint32_t> get_count_cigarotypes();
std::vector<uint32_t> get_nlengths();
std::vector<uint32_t> get_start_positions();
std::vector<std::vector<std::string>> get_qnames();
void transform_mutations();
void set_positions(position pos);
std::vector<position> get_positions();
void populate_positions();
#endif
