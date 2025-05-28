#include <iostream>
#include "primer_bed.h"
#include "genomic_position.h"
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

  ITNode(Interval value)
      : data(new Interval(value)),
        left(nullptr),
        right(nullptr),
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
  void print_amplicons(ITNode *root);
  void get_max_pos(ITNode *root);
  int unpaired_primers(ITNode *root, primer prim);
  void find_read_amplicon(ITNode *root, uint32_t lower, uint32_t upper, ITNode*&node, uint32_t &amp_dist);
  void calculate_overlaps(ITNode *root, std::vector<genomic_position> &positions);

 public:
  IntervalTree();
  uint32_t max_pos=0;
  std::vector<std::vector<uint32_t>> overlaps;
  std::vector<uint32_t> test_test;
  std::vector<uint32_t> flagged_positions; //positions where freq flux occurs MIGHT NOT NEED

  void insert(Interval data) { _root = insert_node_balanced(_root, data); }
  bool is_interval_contained(Interval data) { return is_interval_contained(_root, data); }
  void inOrder() { inOrder(_root); }
  std::string pre_order_with_level() { return pre_order_with_level(_root, 0); }
  void print_amplicons() {print_amplicons(_root);}
  void get_max_pos() {get_max_pos(_root);}
  int unpaired_primers(primer prim) { return unpaired_primers(_root, prim);}
  void calculate_overlaps(std::vector<genomic_position> &positions) {calculate_overlaps(_root, positions);}
  void find_read_amplicon(uint32_t lower, uint32_t upper, ITNode*&node, uint32_t &amp_dist) {find_read_amplicon(_root, lower, upper, node, amp_dist);}
};

int unpaired_primers(ITNode *root, primer prim);
bool node_compare(ITNode *node1, ITNode *node2);
IntervalTree populate_amplicons(std::string pair_info_file, std::vector<primer> &primers);
#endif
