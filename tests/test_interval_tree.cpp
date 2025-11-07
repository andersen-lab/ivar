#include <iostream>

#include "../src/interval_tree.h"
#include "../src/primer_bed.h"

int test_itree_overlap(IntervalTree tree, Interval queries[], int num_tests,
                       bool expected[]) {
  int result = 0;
  for (int i = 0; i < num_tests; i++) {
    result = tree.is_interval_contained(queries[i]);
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


int test_left_rotate() {
  IntervalTree tree = IntervalTree();

  // Right unbalanced tree without rotations
  tree.insert(Interval(10, 15));
  tree.insert(Interval(20, 25));
  tree.insert(Interval(30, 35));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[(20,25)(10,35)](0), [(10,15)(10,15)](1), [(30,35)(30,35)](1), ") ? 0 : 1;
}

int test_right_rotate() {
  IntervalTree tree = IntervalTree();

  // Left unbalanced tree without rotations
  tree.insert(Interval(30, 35));
  tree.insert(Interval(20, 25));
  tree.insert(Interval(10, 15));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[(20,25)(10,35)](0), [(10,15)(10,15)](1), [(30,35)(30,35)](1), ") ? 0 : 1;
}

int test_right_left_rotate() {
  IntervalTree tree = IntervalTree();

  // Right unbalanced tree without rotations
  tree.insert(Interval(10, 15));
  tree.insert(Interval(30, 35));
  tree.insert(Interval(20, 25));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[(20,25)(10,35)](0), [(10,15)(10,15)](1), [(30,35)(30,35)](1), ") ? 0 : 1;
}

int test_left_right_rotate() {
  IntervalTree tree = IntervalTree();

  // Right unbalanced tree without rotations
  tree.insert(Interval(30, 35));
  tree.insert(Interval(10, 15));
  tree.insert(Interval(20, 25));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[(20,25)(10,35)](0), [(10,15)(10,15)](1), [(30,35)(30,35)](1), ") ? 0 : 1;
}

int test_interval_overlap(){
  IntervalTree tree = IntervalTree();

  // Insert several test intervals
  tree.insert(Interval(5, 20));
  tree.insert(Interval(10, 30));
  tree.insert(Interval(15, 25));
  tree.insert(Interval(25, 40));
  tree.insert(Interval(50, 60));

  std::string pre_order_str = tree.pre_order_with_level();

  if(pre_order_str != "[(10,30)(5,60)](0), [(5,20)(5,20)](1), [(25,40)(15,60)](1), [(15,25)(15,25)](2), [(50,60)(50,60)](2), ")
    return 1;


  Interval ints[] = {Interval(15, 25), Interval(16, 24), Interval(5, 60), Interval(31, 35), Interval(45, 55), Interval(70, 80)};
  bool results[] = {true, true, false, true, false, false};

  for(int i = 0; i < 5; i++){
    if(results[i] != tree.is_interval_contained(ints[i]))
      return 1;
  }

  // Find if read overlaps fully with any interval
  std::vector<std::vector<int>> test_intervals = {
      {12, 18},
      {6, 19},
      {16, 27},
      {55, 59},
      {17, 19}
  };

  // expected nodes output should be in order of traversal of interval tree
  std::vector<std::vector<Interval>> expected_nodes = {
      {Interval(10, 30), Interval(5, 20)},
      {Interval(5, 20)},
      {Interval(10, 30)},
      {Interval(50, 60)},
      {Interval(10, 30), Interval(5, 20), Interval(15, 25)}
  };

  for(int i = 0; i < test_intervals.size(); i++){
    std::vector<ITNode*> nodes;
    tree.find_read_amplicon(test_intervals[i][0], test_intervals[i][1], nodes);
    if (nodes.size() != expected_nodes[i].size()) {
      std::cerr << "Nodes size did not match for test interval " << test_intervals[i][0] << "-" << test_intervals[i][1] << "\n";
      return 1;
    }
    for(int j = 0; j < nodes.size(); j++){
        if(nodes[j]->data->low != expected_nodes[i][j].low || nodes[j]->data->high != expected_nodes[i][j].high) {
          std::cerr << "Node interval did not match for test interval " << test_intervals[i][0] << "-" << test_intervals[i][1] << "\n";
          return 1;
        }
    }
  }

  return 0;
}

int main() {
  int result = 0;

  result += test_left_rotate();
  result += test_right_rotate();
  result += test_left_right_rotate();
  result += test_right_left_rotate();
  result += test_interval_overlap();

  return result;
}