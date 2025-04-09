#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <sstream>
#include <algorithm>

#ifndef parse_gff
#define parse_gff

/*
GFF3 file format:
Defined at
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

| Col        | Desc                                                          |
| seqid      | ID of ref                                                     |
| source     | algorithm/operating procedure                                 |
| type       | type of feature                                               |
| start      | 1-based start                                                 |
| end        | start <= end                                                  |
| score      | score of feature                                              |
| strand     | +/-                                                           |
| phase      | Number of bases to be removed with reference to reading frame |
| attributes | list of attributes in the format tag=value;tag-value..        |
*/

/*
Custom GFF3 attributes for RNA editing

EditPosition: Position to make insertion
EditSequence: Sequence to be inserted

*/

const std::string EDIT_POSITION = "EditPosition";
const std::string EDIT_SEQUENCE = "EditSequence";

class gff3_feature {
 public:
  gff3_feature(std::string line);
  int print();

  std::string get_seqid();
  std::string get_source();
  std::string get_type();
  uint64_t get_start();
  uint64_t get_end();
  char get_strand();
  int get_phase();
  std::map<std::string, std::string> get_attributes();
  std::string get_attribute(std::string key);
  gff3_feature* get_previous();
  gff3_feature* get_next();

  int set_seqid();
  int set_source();
  int set_type();
  int set_start();
  int set_end();
  int set_strand();
  int set_phase();
  int set_attributes(std::string attr);
  int64_t get_edit_position();
  std::string get_edit_sequence();
  void set_previous(gff3_feature* prev);
  void set_next(gff3_feature* next);

 private:
  std::string seqid, source, type;
  std::map<std::string, std::string> attributes;
  uint64_t start, end;
  float score;
  char strand;
  int phase;
  gff3_feature* previous_feature;
  gff3_feature* next_feature;
};

class gff3 {
 public:
  gff3();
  gff3(std::string path);
  std::vector<gff3_feature> get_features();
  int print();
  int read_file(std::string path);
  std::vector<gff3_feature> query_features(uint64_t pos, std::string type);
  int get_count();
  bool empty();
  void link_features();

 private:
  std::vector<gff3_feature> features;
  // Flag to see if file has been populated
  bool is_empty;
};

#endif
