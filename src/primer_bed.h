#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>

#ifndef primer_bed
#define primer_bed

class primer {
 private:
  std::string region;
  uint32_t start;		// 0 based
  uint32_t end;			// 0 based
  std::string name;
  int score;
  char strand;
  int16_t pair_indice;
  int16_t indice;
  uint32_t read_count;

 public:
  std::string get_name();
  std::string get_region();
  int get_score();
  uint32_t get_start() const;
  uint32_t get_end() const;
  char get_strand();
  int get_length();
  int16_t get_pair_indice();
  int16_t get_indice() const;
  uint32_t get_read_count() const;
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
  bool operator == (const primer& p) const {
    return (indice == p.get_indice()) ? true : false;
  }

};

std::vector<primer> populate_from_file(std::string path, int32_t offset);
std::vector<primer> populate_from_file(std::string path);
std::vector<primer> get_primers(std::vector<primer> p, unsigned int pos);
int get_primer_indice(std::vector<primer> p, std::string name);
int populate_pair_indices(std::vector<primer> &primers, std::string path);
primer get_min_start(std::vector<primer> primers);
primer get_max_end(std::vector<primer> primers);

#endif
