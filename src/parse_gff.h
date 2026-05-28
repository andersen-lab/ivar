#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <sstream>
#include <vector>

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
  std::string get_type() const;
  uint64_t get_start() const;
  uint64_t get_end() const;
  char get_strand() const;
  int get_phase() const;
  std::map<std::string, std::string> get_attributes();
  std::string get_attribute(std::string key) const;

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

 private:
  std::string seqid, source, type;
  std::map<std::string, std::string> attributes;
  uint64_t start, end;
  float score;
  char strand;
  int phase;
};

// Combine CDSs with same ID into a group
class cds_group {
 public:
  struct codon_triple {
    int64_t g0, g1, g2;
    int64_t min_pos, max_pos;
  };

  cds_group();
  void add_segment(const gff3_feature &f);
  void sort_and_finalize_segments();

  // pos is 1-based genomic; cds_pos is 0-based offset in 5' to 3'
  // (for - strand, cds_pos=0 corresponds to the highest genomic coord).
  int64_t genomic_to_cds_pos(int64_t pos) const; // -1 if not in group
  int64_t cds_to_genomic_pos(int64_t cds_pos) const; // -1 if past length

  char get_strand() const;
  int get_phase() const; // phase of 5'-most segment
  int64_t length() const; // sum of segment lengths
  bool contains(int64_t pos) const;
  const std::string &get_id() const;
  const std::string &get_gene() const;
  const std::vector<gff3_feature> &segments() const;
  const std::vector<codon_triple> &codon_triples() const;

 private:
  std::vector<gff3_feature> segments_;
  std::vector<int64_t> cumulative_len_before_;
  std::vector<int64_t> cds_genomic_index_;
  std::vector<codon_triple> codon_triples_;
  std::string id_;
  std::string gene_;
  char strand_;
};

class gff3 {
 public:
  gff3();
  gff3(std::string path);
  std::vector<gff3_feature> get_features();
  int print();
  int read_file(std::string path);
  std::vector<gff3_feature> query_features(uint64_t pos, std::string type);
  std::vector<cds_group> query_cds_groups(uint64_t pos) const;
  const std::vector<cds_group> &get_cds_groups() const;
  const cds_group *find_cds_group_by_id(const std::string &id) const;
  int get_count();
  bool empty();

 private:
  void build_cds_groups();
  std::vector<gff3_feature> features;
  std::vector<cds_group> cds_groups_;
  std::map<std::string, size_t> cds_group_index_;
  // Flag to see if file has been populated
  bool is_empty;
};

#endif
