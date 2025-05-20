#include <iostream>

#include "primer_bed.h"
using namespace std;

#ifndef interval_tree
#define interval_tree

// Interval
class Interval{
 public:
  int low, high;

  Interval(int val1, int val2)
      : low(std::min(val1, val2)), high(std::max(val1, val2)) {}  // constructor
};

// A node in IntervalTree
class ITNode {
 public:
  void set_haplotypes(primer prim);
  Interval *data; // pointer to node's interval data object
  ITNode *left, *right; // pointer to node's left & right child node objects
  int max;
  int height; // height of node
  std::vector<position> amp_positions;  //data for every position on amplicon

  ITNode(Interval value)
      : data(new Interval(value)),
        left(nullptr),
        right(nullptr),
 improved_interval_tree
        max(value.high),
        height(0){}

  // Get method to handle null case
  int get_node_height(ITNode* node) const {
    return (node != nullptr) ? node->height : -1;
  }

  void update_height() {
    height = 1 + std::max(get_node_height(left), get_node_height(right));
  }

  int get_balance() const {
    return get_node_height(left) - get_node_height(right);
  }

  void update_max() {
    max = data->high;
    if (left != nullptr)
      max = std::max(max, left->max);
    if (right != nullptr)
      max = std::max(max, right->max);
  }
};

/////////////////////////////////////////////////////////////////////////////////////////
// IntervalTree class
class IntervalTree {
 private:
  ITNode *_root;
  ITNode *left_rotate(ITNode *x);
  ITNode* right_rotate(ITNode *y);
  ITNode* insert_node_balanced(ITNode *node, Interval data);

  void insert(ITNode *root, Interval data);
  bool is_interval_contained(ITNode *root, Interval data);
  void inOrder(ITNode *root);
  std::string pre_order_with_level(ITNode *root, int level);
  void amplicon_position_pop(ITNode *root);
  void print_amplicons(ITNode *root);
  void get_max_pos(ITNode *root);
  int unpaired_primers(ITNode *root, primer prim);
  void combine_haplotypes(ITNode *root, uint32_t &counter);
  void write_out_frequencies(ITNode *root, std::string filename);
  void detect_abberations(ITNode *root, uint32_t pos);
  void find_read_amplicon(ITNode *root, uint32_t lower, uint32_t upper, bool &found, std::string read_name, uint32_t &amp_start, uint32_t &amp_dist);
  void assign_read_amplicon(ITNode *root, uint32_t amp_start, std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities, uint8_t min_qual);

 public:
  IntervalTree();
  uint32_t max_pos=0;
  std::vector<std::vector<uint32_t>> overlaps;
  std::vector<position> test_flux; //storage for looking at pos across all amps
  std::vector<uint32_t> test_test;
  std::unordered_map<uint32_t, position> variants; //all variants across every position                                 
  std::unordered_map<uint32_t, position> amp_positions;
  std::vector<uint32_t> flagged_positions; //positions where freq flux occurs MIGHT NOT NEED

  void insert(Interval data) { _root = insert_node_balanced(_root, data); }
  bool is_interval_contained(Interval data) { return is_interval_contained(_root, data); }
  void inOrder() { inOrder(_root); }
  std::string pre_order_with_level() { return pre_order_with_level(_root, 0); }
  void print_amplicons() {print_amplicons(_root);}
  void get_max_pos() {get_max_pos(_root);}
  int unpaired_primers(primer prim) { return unpaired_primers(_root, prim);}
  void detect_abberations(uint32_t pos) {detect_abberations(_root, pos);}
  void combine_haplotypes(uint32_t &counter) {combine_haplotypes(_root, counter);}
  void write_out_frequencies(std::string filename){write_out_frequencies(_root, filename);}
  void populate_variants(uint32_t last_position);
  void add_read_variants(std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities, uint8_t min_qual);
  void find_read_amplicon(uint32_t lower, uint32_t upper, bool &found, std::string read_name, uint32_t &amp_start, uint32_t &amp_dist) {find_read_amplicon(_root, lower, upper, found, read_name, amp_start, amp_dist);}
  void assign_read_amplicon(uint32_t amp_start, std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities, uint8_t min_qual) {assign_read_amplicon(_root, amp_start, positions, bases, qualities, min_qual);}
  void amplicon_position_pop() {amplicon_position_pop(_root);}
};

int unpaired_primers(ITNode *root, primer prim);
IntervalTree populate_amplicons(std::string pair_info_file, std::vector<primer> &primers);
IntervalTree amplicon_position_pop();
#endif
