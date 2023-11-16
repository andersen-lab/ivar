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
  std::vector<position> amp_positions;

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
  void get_size(ITNode *root);
  void set_haplotypes(ITNode *root, primer prim);
  public:
  uint32_t count;
  IntervalTree();  // constructor
  void insert(Interval data) { insert(_root, data); }
  bool envelopSearch(Interval data) { return envelopSearch(_root, data); }
  void inOrder() { inOrder(_root); }
  void print_amplicons() {print_amplicons(_root);}
  void get_size() {get_size(_root);}
  void set_haplotypes (primer prim) {set_haplotypes(_root, prim);}
};

void get_size();
void set_haplotypes(ITNode *root, primer prim);
IntervalTree populate_amplicons(std::string pair_info_file,
                                std::vector<primer> &primers);
#endif
