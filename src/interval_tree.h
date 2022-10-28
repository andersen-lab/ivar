#include <iostream>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "primer_bed.h"
#include "allele_functions.h"
using namespace std;

#ifndef interval_tree
#define interval_tree

// Structure to represent an interval
class Interval{   public:
  Interval(uint32_t val1, uint32_t val2, uint32_t val3, uint32_t val4): low(std::min(val1, val2)), high(std::max(val1, val2)), low_inner(val3), high_inner(val4) {}  // constructor
  uint32_t low, high, low_inner, high_inner;
};
// Structure to represent a node in Interval Search Tree
class ITNode{
  /*
    public:
    ITNode(Interval *values): data(value), left(nullptr), right(nullptr) {}  // constructor
    int max;
    // Getters - access member functions
    Interval getData()const;
    ITNode getLeft()const;
    ITNode getRight()const;
    // Setters - access member functions
    void setLeft(ITNode *node);
    void setRight(ITNode *node);
  */
public:
  ITNode(Interval value): data(new Interval(value)), left(nullptr), right(nullptr), max(value.high) {}  // constructor
  Interval *data;  // pointer to node's interval data object
  ITNode *left, *right; // pointer to node's left & right child node objects
  uint32_t max;
  std::vector<std::vector<int>> haplotypes;
  std::vector<std::vector<uint32_t>> positions;
  std::vector<std::vector<uint32_t>> ranges; //the starts and ends of the reads

  std::vector<std::vector<int>> final_haplotypes; //these are summary for the amplicon
  std::vector<uint32_t> final_positions;
  std::vector<float> frequency;

  std::vector<uint32_t> suspect_positions;
  std::vector<uint32_t> suspect_haplotypes;

  double reverse=0; //these are added to track how many reverse/forward reads have mutations
  double forward=0;
  double mut_reverse=0;
  double mut_forward=0;
  int read_count=0;

};

/////////////////////////////////////////////////////////////////////////////////////////
// IntervalTree class
class IntervalTree{
private:
  ITNode *_root;
  void insert(ITNode *root, Interval data);
  bool envelopSearch(ITNode *root, Interval data);
  void inOrder(ITNode * root); //prints amplicons starts and ends
  void print_amplicon_summary(ITNode *root); //prints summary of unique haplotypes and frequencies
  void clear(ITNode *root); //removes all information from each node
  void dump_amplicon_summary(ITNode *root, std::string filename); //dump amplicon summaries to json file
  void find_amplicon_per_read(ITNode *root, uint32_t start, uint32_t end, std::vector<int> haplotypes, 
      std::vector<uint32_t> positions, bool reverse, std::vector<uint32_t> ranges, std::vector<position> &all_positions, bool primer_mutation, bam1_t *r);

public:
  IntervalTree();  // constructor
  void insert(Interval data){ insert(_root, data);}
  bool envelopSearch(Interval data){ return envelopSearch(_root, data);}
  void inOrder() {inOrder(_root);}
  void print_amplicon_summary() {print_amplicon_summary(_root);}
  void clear() {clear(_root);}
  void dump_amplicon_summary(std::string filename) {dump_amplicon_summary(_root, filename);}
  ITNode *iterate_nodes(ITNode *root); //used to returns nodes iteratively
  ITNode *iterate_nodes(){return iterate_nodes(_root);}
  void find_amplicon_per_read(uint32_t start, uint32_t end, std::vector<int> haplotypes, 
      std::vector<uint32_t> positions, bool reverse, std::vector<uint32_t> ranges, std::vector<position> &all_positions, bool primer_mutation, bam1_t *r){find_amplicon_per_read(_root, start, end, haplotypes, positions, reverse, ranges, all_positions, primer_mutation, r);}
};

IntervalTree populate_amplicons(std::string pair_info_file, std::vector<primer> &primers);
void remove_low_quality_nts(ITNode *node, std::vector<position> all_positions);
#endif
