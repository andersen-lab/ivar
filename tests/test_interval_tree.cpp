#include <iostream>

#include "../src/interval_tree.h"
#include "../src/primer_bed.h"

int test_itree_overlap(IntervalTree tree, Interval queries[], int num_tests,
                       bool expected[]) {
  int result = 0;
  for (int i = 0; i < num_tests; i++) {
    result = tree.envelopSearch(queries[i]);
    if (result != expected[i]) {
      std::cout << "Interval Tree overlap behavior incorrect for interval "
                << queries[i].low << ":" << queries[i].high << " - "
                << "Expected: " << expected[i] << "Got: " << result
                << std::endl;
      return 1;
    }
  }
  return 0;
}

int main() {
  int result = 0;
  Interval ints[6] = {Interval(15, 20), Interval(30, 10), Interval(17, 19),
                      Interval(5, 20),  Interval(12, 15), Interval(30, 40)};
  int n = sizeof(ints) / sizeof(ints[0]);
  IntervalTree tree = IntervalTree();
  // populate interval tree
  for (int i = 0; i < n; i++) {
    tree.insert(ints[i]);
  }
  std::cout << "In order traversal of Interval Tree:" << std::endl;
  tree.inOrder();

  Interval queries[4] = {Interval(15, 20), Interval(9, 30), Interval(31, 38),
                         Interval(7, 22)};
  bool expected[4] = {true, false, true, false};
  int num_tests = sizeof(queries) / sizeof(queries[0]);
  result = test_itree_overlap(tree, queries, num_tests, expected);
  return result;
}