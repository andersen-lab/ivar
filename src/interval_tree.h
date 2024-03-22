#include <iostream>

#include "primer_bed.h"
using namespace std;

#ifndef interval_tree
#define interval_tree

// Structure to represent an interval
class Interval{
 public:
  Interval(int val1, int val2)
      : low(std::min(val1, val2)), high(std::max(val1, val2)) {}  // constructor
  int low, high;
};

// Structure to represent a node in Interval Search Tree
class ITNode {
  /*
    public:
    ITNode(Interval *values): data(value), left(nullptr), right(nullptr) {}  //
    constructor int max;
    // Getters - access member functions
    Interval getData()const;
    ITNode getLeft()const;
    ITNode getRight()const;
    // Setters - access member functions
    void setLeft(ITNode *node);
    void setRight(ITNode *node);
  */
  public:
  void set_haplotypes(primer prim);
  ITNode(Interval value)
      : data(new Interval(value)),
        left(nullptr),
        right(nullptr),
        max(value.high) {}  // constructor
  Interval *data;           // pointer to node's interval data object
  ITNode *left, *right;     // pointer to node's left & right child node objects
  int max;
  std::vector<position> amp_positions;  //data for every position on amplicon                             
};

/////////////////////////////////////////////////////////////////////////////////////////
// IntervalTree class
class IntervalTree {
 private:
  ITNode *_root;
  void insert(ITNode *root, Interval data);
  bool envelopSearch(ITNode *root, Interval data);
  void inOrder(ITNode *root);
  void print_amplicons(ITNode *root);
  void get_max_pos(ITNode *root);
  void set_haplotypes(ITNode *root, primer prim);
  int unpaired_primers(ITNode *root, primer prim);
  void combine_haplotypes(ITNode *root);
  void detect_abberations(ITNode *root, uint32_t pos);
  void detect_amplicon_overlaps(ITNode *root, uint32_t pos);
  void detect_primer_issues(ITNode *root, uint32_t pos);
  public:
  uint32_t max_pos=0;
  std::vector<std::vector<uint32_t>> overlaps;
  std::vector<position> test_flux; //storage for looking at pos across all amps
  std::vector<position> variants; //all variants across every position                                 
  std::vector<uint32_t> flagged_positions; //positions where freq flux occurs MIGHT NOT NEED
  IntervalTree();  // constructor
  void insert(Interval data) { insert(_root, data); }
  bool envelopSearch(Interval data) { return envelopSearch(_root, data); }
  void inOrder() { inOrder(_root); }
  void print_amplicons() {print_amplicons(_root);}
  void get_max_pos() {get_max_pos(_root);}
  void set_haplotypes(primer prim) {set_haplotypes(_root, prim);}
  int unpaired_primers(primer prim) { return unpaired_primers(_root, prim);}
  void detect_abberations(uint32_t pos) {detect_abberations(_root, pos);}
  void detect_amplicon_overlaps(uint32_t pos) {detect_amplicon_overlaps(_root, pos);}
  void detect_primer_issues(uint32_t pos) {detect_primer_issues(_root, pos);}
  void combine_haplotypes() {combine_haplotypes(_root);}
  void add_read_variants(uint32_t *cigar, uint32_t start_pos, uint32_t nlength, uint8_t *sequence, uint8_t *aux, uint8_t *quality, std::string qname);
  void populate_variants();
  
};

void combine_haplotypes();
void detect_abberations(ITNode *root, uint32_t find_position);
void get_max_pos();
void set_haplotypes(ITNode *root, primer prim);
void add_read_variants(uint32_t *cigar, uint32_t start_pos, uint32_t nlength, uint8_t *sequence, uint8_t *aux, uint8_t *quality, std::string qname);
void populate_variants();
int unpaired_primers(ITNode *root, primer prim);
void detect_primer_issues(ITNode *root, uint32_t find_position);
void detect_amplicon_overlaps(ITNode *root, uint32_t find_position);
IntervalTree populate_amplicons(std::string pair_info_file,
                                std::vector<primer> &primers);
#endif
