#include "interval_tree.h"
#include <unordered_map>

// Constructor for initializing an Interval Tree
IntervalTree::IntervalTree() { _root = NULL; }
std::vector<allele> add_allele_vectors(std::vector<allele> new_alleles, std::vector<allele> return_alleles){
  /*
   * @param return_alleles : the alleles we're saving to the amplicon
   * @param new_alleles : the alleles from the primer
   */
  for(uint32_t i=0; i < new_alleles.size(); i++){
    bool found = false;
    for(uint32_t j=0; j < return_alleles.size(); j++){
      if ((return_alleles[j].nuc == new_alleles[i].nuc) && new_alleles[i].depth > 0){
        return_alleles[j].depth += new_alleles[i].depth;
        return_alleles[j].mean_qual += new_alleles[i].mean_qual;
        found = true;
        break;
      }
    }
    //we don't have this allele in our amplicon haplotype
    if (!found) {
      return_alleles.push_back(new_alleles[i]);
    }
  }
  return(return_alleles);
}

void IntervalTree::write_out_frequencies(ITNode *root, std::string filename){
  //POS\tALLELE\tDEPTH\tFREQ\tAMP_START\tAMP_END\n
  if (root==NULL) return;
  std::fstream file(filename, std::ios::in | std::ios::out | std::ios::app);
  for(uint32_t i=0; i < root->amp_positions.size(); i++){
    if(root->amp_positions[i].depth == (uint32_t)0) continue;
    uint32_t total_depth = 0;
    for(uint32_t j=0; j < root->amp_positions[i].alleles.size(); j++){
      //if(root->amp_positions[i].alleles[j].nuc != "-"){
        total_depth += root->amp_positions[i].alleles[j].depth;
      //}
    }
    for(uint32_t j=0; j < root->amp_positions[i].alleles.size(); j++){
      if(root->amp_positions[i].alleles[j].depth == 0) continue;
      file << std::to_string(root->amp_positions[i].pos) << "\t";
      file << root->amp_positions[i].alleles[j].nuc << "\t";
      file << std::to_string(root->amp_positions[i].alleles[j].depth) << "\t";
      file << std::to_string((double)root->amp_positions[i].alleles[j].depth/(double)total_depth) << "\t";
      file << std::to_string(root->data->low) << "\t";
      file << std::to_string(root->data->high) << "\n";

    }
  } 
  file.close();
  write_out_frequencies(root->right, filename);
}

void IntervalTree::combine_haplotypes(ITNode *root, uint32_t &counter){
  if (root == NULL) return;
  //traverse left subtree
  combine_haplotypes(root->left, counter);

  //process current node
  for (auto [pos_key, amp_pos] : root->amp_positions) {
    if (amp_pos.depth == 0) continue;

    //update position depth
    position& var = variants[pos_key];
    var.depth += amp_pos.depth;

    //update alleles
    var.alleles = add_allele_vectors(amp_pos.alleles, var.alleles);
    var.amplicon_numbers.push_back(counter);

    double inv_depth = 1.0 / amp_pos.depth;
    for (const auto& allele : amp_pos.alleles) {
      var.amplicon_frequencies[allele.nuc].push_back(allele.depth * inv_depth);
    }
  }

  ++counter;

  //traverse right subtree
  combine_haplotypes(root->right, counter);
}


void IntervalTree::assign_read_amplicon(ITNode *node, std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities, uint8_t min_qual) {
  for (uint32_t i = 0; i < positions.size(); i++) {
    if (qualities[i] < (uint32_t)min_qual) continue;
      node->amp_positions[positions[i]].update_alleles(bases[i], 1, qualities[i]);
  }
}

void IntervalTree::find_read_amplicon(ITNode *root, uint32_t lower, uint32_t upper, ITNode* &node, uint32_t &amp_dist) {
  if (root == NULL) return;

  //check if current node's interval fully contains [lower, upper]
  if ((uint32_t)root->data->low <= lower && upper <= (uint32_t)root->data->high) {
    uint32_t dist = (lower - root->data->low) + (root->data->high - upper);
    if (dist < amp_dist) {
      amp_dist = dist;
      node = root;
    }
  }

  //traverse left if there's any chance of finding a containing interval
  if (root->left && lower <= (uint32_t)root->left->data->high) {
    find_read_amplicon(root->left, lower, upper, node, amp_dist);
  }

  //traverse right if there's any chance of finding a containing interval
  if (root->right && upper >= (uint32_t)root->right->data->low) {
    find_read_amplicon(root->right, lower, upper, node, amp_dist);
  }
}


void IntervalTree::amplicon_position_pop(ITNode *root) {
  if (root == NULL) return;

  // Traverse left subtree
  amplicon_position_pop(root->left);

  // Populate amp_positions for current node
  for (uint32_t i = root->data->low; i < (uint32_t)root->data->high; i++) {
    position add_pos;
    add_pos.pos = i;
    add_pos.alleles = populate_basic_alleles();
    root->amp_positions[i] = add_pos;
  }

  // Traverse right subtree
  amplicon_position_pop(root->right);
}

void IntervalTree::detect_abberations(ITNode *root, uint32_t find_position) {
  if (root == NULL) return;

  //check current node for a matching position
  auto it = root->amp_positions.find(find_position);
  if (it != root->amp_positions.end() && it->second.depth > 0) {
    position& amp_pos = it->second;
    test_flux.push_back(amp_pos);
    test_test.push_back(root->data->low); //track which amplicon/interval matched
  }

  //traverse both sides of the tree to catch all overlaps
  detect_abberations(root->left, find_position);
  detect_abberations(root->right, find_position);
}


void IntervalTree::add_read_variants(std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities, uint8_t min_qual){
  for(uint32_t i=0; i < positions.size(); i++){
    if(qualities[i] >= (uint32_t)min_qual){
      variants[positions[i]].update_alleles(bases[i], 1, qualities[i]);  
    }
  }
}

void IntervalTree::populate_variants(uint32_t last_position){
  for(uint32_t i=0; i <= last_position; i++){
    position tmp;
    tmp.alleles = populate_basic_alleles();
    tmp.pos = i;
    variants[i] = tmp;
  }
}

int IntervalTree::unpaired_primers(ITNode *root, primer prim){
  if (root==NULL) return 0;
  char strand = prim.get_strand();
  if (strand == '+' && ((int)prim.get_start() == root->data->low)){
    return 1;
  } else if (strand == '-' && ((int)prim.get_end()+1 == root->data->high)){
    return 1;
  }
  int ret = unpaired_primers(root->right, prim);
  return ret;
}

// [Deprecated] Use insert_node_balanced()
// A utility function to insert a new Interval Search Tree Node
// This is similar to BST Insert.  Here the low value of interval
// is used tomaintain BST property
void IntervalTree::insert(ITNode *root, Interval data) {
  // Base case: Tree is empty, new node becomes root
  if (root == NULL) {
    root = new ITNode(data);
    _root = root;
  } else {
    // Get low value of interval at root
    int l = root->data->low;
    // If root's low value is greater, then new interval goes to
    // left subtree
    if (data.low < l) {
      if (!root->left) {
        ITNode *tmpNode = new ITNode(data);
        // std::cout << data.low << ":" << data.high << "->insertLeft" <<
        // std::endl;
        root->left = tmpNode;
      } else {
        insert(root->left, data);
      }
    } else {
      if (!root->right) {
        ITNode *tmpNode = new ITNode(data);
        // std::cout << data.low << ":" << data.high << "->insertRight" <<
        // std::endl;
        root->right = tmpNode;
      } else {
        insert(root->right, data);
      }
    }
  }
  // update max value of ancestor node
  if (root->max < data.high) root->max = data.high;
}

// A utility function to check if the 1st interval envelops the second
bool check_interval_contained(Interval i1, Interval i2) {
  if (i1.low <= i2.low && i1.high >= i2.high) return true;
  return false;
}

// The main function that searches an interval i in a given
// Interval Tree.
bool IntervalTree::is_interval_contained(ITNode *root, Interval i) {
  // Base Case, tree is empty
  // std::cout << root->data->low << ":" << root->data->high << std::endl;
  if (root == NULL) return false;

  // If given interval overlaps with root
  if (check_interval_contained(*(root->data), i)) return true;

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // be enveloped by an amplicon in left subtree
  if (root->left != NULL && root->left->max >= i.high)
    return is_interval_contained(root->left, i);

  // Else interval can only be enveloped by amplicon in right subtree
  return is_interval_contained(root->right, i);
}

void IntervalTree::get_max_pos(ITNode *root){
  if (root == NULL) return;
  if (root->data->high > (int)max_pos) {
    max_pos = (uint32_t) root->data->high;
  }
  get_max_pos(root->right);
}

// A helper function for inorder traversal of the tree
void IntervalTree::inOrder(ITNode *root) {
  if (root == NULL) return;
  inOrder(root->left);
  std::cerr << "[" << root->data->low << ", " << root->data->high << "]"
            << " max = " << root->max << endl;
  inOrder(root->right);
}

std::string IntervalTree::pre_order_with_level(ITNode *root,  int level) {
  if (root == nullptr) return "";
  std::string pre_order_str;
  pre_order_str = "[" + std::to_string(root->data->low) + "," + std::to_string(root->data->high) + "]";
  pre_order_str += "(" + std::to_string(level) + "), ";
  pre_order_str += pre_order_with_level(root->left, level + 1);
  pre_order_str += pre_order_with_level(root->right, level + 1);
  return pre_order_str;
}

// A stand-alone function to create a tree containing the coordinates of each
// amplicon based on user-specified primer pairs
IntervalTree populate_amplicons(std::string pair_info_file,
                                std::vector<primer> &primers) {
  int amplicon_start = -1;
  int amplicon_end = -1;
  IntervalTree tree = IntervalTree();
  populate_pair_indices(primers, pair_info_file);
  for (auto &p : primers) {
    if (p.get_strand() == '+') {
      if (p.get_pair_indice() != -1) {
        amplicon_start = p.get_start();
        amplicon_end = primers[p.get_pair_indice()].get_end() + 1;
        tree.insert(Interval(amplicon_start, amplicon_end));
      }
    }
  }
  return tree;
}

void IntervalTree::print_amplicons(ITNode *root){
  if (root==NULL) return;
  std::cout << "\nLow:" << root->data->low << " High: " << root->data->high << std::endl;
  print_amplicons(root->right);
}

// node is right unbalanced
ITNode* IntervalTree::left_rotate(ITNode *node) {
  ITNode *right = node->right;
  ITNode *right_left = right->left;

  // Rotate
  right->left = node;
  node->right = right_left;

  node->update_height();
  right->update_height();
  node->update_max();
  right->update_max();

  // New root of subtree
  return right;
}

// node if left unbalanced
ITNode* IntervalTree::right_rotate(ITNode *node) {
  ITNode *left = node->left;
  ITNode *left_right = left->right;

  // Rotate
  left->right = node;
  node->left = left_right;

  node->update_height();
  left->update_height();
  node->update_max();
  left->update_max();

  // New root of subtree
  return left;
}

ITNode* IntervalTree::insert_node_balanced(ITNode *node, Interval data) {
  if (node == nullptr)
    return new ITNode(data);

  if (data.low < node->data->low)
    node->left = insert_node_balanced(node->left, data); // Insert into left subtree
  else
    node->right = insert_node_balanced(node->right, data); // Insert into right subtree

  node->update_height();
  node->update_max();

  int balance = node->get_balance();

  if(balance > 1) { // Left unbalanced
    if(data.low < node->left->data->low) // Node inserted in left left
      return right_rotate(node);
    else { // Node inserted in left right
      node->left = left_rotate(node->left);
      return right_rotate(node);
    }
  }

  if(balance < -1) {
    if(data.low >= node->right->data->low) // Node inserted in right right
      return left_rotate(node);
    else { // Node inserted in right left
      node->right = right_rotate(node->right);
      return left_rotate(node);
    }
  }

  return node;
}
