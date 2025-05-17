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

  return (pre_order_str == "[20,25](0), [10,15](1), [30,35](1), ") ? 0 : 1;
}

int test_right_rotate() {
  IntervalTree tree = IntervalTree();

  // Left unbalanced tree without rotations
  tree.insert(Interval(30, 35));
  tree.insert(Interval(20, 25));
  tree.insert(Interval(10, 15));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[20,25](0), [10,15](1), [30,35](1), ") ? 0 : 1;
}

int test_right_left_rotate() {
  IntervalTree tree = IntervalTree();

  // Right unbalanced tree without rotations
  tree.insert(Interval(10, 15));
  tree.insert(Interval(30, 35));
  tree.insert(Interval(20, 25));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[20,25](0), [10,15](1), [30,35](1), ") ? 0 : 1;
}

int test_left_right_rotate() {
  IntervalTree tree = IntervalTree();

  // Right unbalanced tree without rotations
  tree.insert(Interval(30, 35));
  tree.insert(Interval(10, 15));
  tree.insert(Interval(20, 25));

  std::string pre_order_str = tree.pre_order_with_level();

  return (pre_order_str == "[20,25](0), [10,15](1), [30,35](1), ") ? 0 : 1;
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

  if(pre_order_str != "[10,30](0), [5,20](1), [25,40](1), [15,25](2), [50,60](2), ")
    return 1;

  Interval ints[] = {Interval(15, 25), Interval(16, 24), Interval(5, 60), Interval(31, 35), Interval(45, 55), Interval(70, 80)};
  bool results[] = {true, true, false, true, false, false};

  for(int i = 0; i < 5; i++){
    if(results[i] != tree.is_interval_contained(ints[i]))
      return 1;
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