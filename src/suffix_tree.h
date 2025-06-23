#include "algorithm"
#include "iostream"
#include "vector"

#ifndef suffix_tree
#define suffix_tree

const std::string alphabet = "ATGCRYSWKMBDHVN-X#@";
const unsigned int MAX_CHAR = alphabet.length();
const int MAX_SIZE = 300;
const int MIN_LENGTH = 3;
const int MAX_MISMATCHES = 2;

class suffix_node {
 public:
  int begin, nchildren, *end;
  suffix_node **children, *parent, *suffix_link;
  suffix_node(int b, int *e, suffix_node *p, suffix_node *l);
  bool is_leaf_node();
  bool contains_child(int ext);
  int get_length();
  int get_depth();
  std::string get_longest_common_substring(std::string s1, std::string s2);
  std::string get_path(std::string s);
  void extend_path(int *e);
  suffix_node *add_child(int ext, int b, int *e, suffix_node *l);
  void add_child(suffix_node *c, int ext);
  suffix_node *get_child(int ext);
  bool contains_depth(int depth);
  void print(std::string s);
  bool walk_next(int &beg, int &suffix_length);
};

suffix_node *build_suffix_tree(std::string s);
std::string get_reverse_complement(std::string rev_read);
std::vector<std::string> read_adapters_from_fasta(std::string p, std::string n);
int trim_adapter(std::string f1, std::string f2, std::string adp_path,
                 std::string p);

#endif
