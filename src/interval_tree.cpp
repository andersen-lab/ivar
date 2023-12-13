#include "interval_tree.h"

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
      if (return_alleles[j].nuc == new_alleles[i].nuc){
        return_alleles[j].depth += new_alleles[i].depth;
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


void IntervalTree::combine_haplotypes(ITNode *root){
  if (root==NULL) return;
  for(uint32_t i=0; i < root->amp_positions.size(); i++){
    uint32_t exists = check_position_exists(root->amp_positions[i].pos, variants);
    //does exist
    if (exists){
      variants[exists].depth += root->amp_positions[i].depth;
      std::vector<allele> new_alleles = add_allele_vectors(root->amp_positions[i].alleles, variants[exists].alleles);
      variants[exists].alleles = new_alleles;
   } else {
      variants.push_back(root->amp_positions[i]);
    }
  }
  combine_haplotypes(root->right);
}

void IntervalTree::detect_abberations(ITNode *root, uint32_t find_position){
  if (root==NULL) return;
  if (find_position < (uint32_t)root->data->low) return;
  if(((uint32_t)root->data->low < find_position) && (find_position < (uint32_t)root->data->high)){
    for(uint32_t i=0; i < root->amp_positions.size(); i++){
      if(find_position == root->amp_positions[i].pos){
        test_flux.push_back(root->amp_positions[i]);
        break;
      }
    }
  }
  detect_abberations(root->right, find_position);

}

void IntervalTree::set_haplotypes(ITNode *root, primer prim){
  if (root==NULL) return;
  char strand = prim.get_strand();
  //these are the bounds on the amplicon
  if(strand == '+' && ((int)prim.get_start() != root->data->low)){
    set_haplotypes(root->right, prim);
  } else if (strand == '-' && ((int)prim.get_end()+1 != root->data->high)){
    set_haplotypes(root->right, prim);
  } else {
    // we found the matching amplion, now we add this cigarotype to the amplicon     
    std::vector<position> tmp_pos = prim.get_positions();
    for (uint32_t i=0; i < tmp_pos.size(); i++){
      position add_pos = tmp_pos[i];
      bool found = false;
      // check if we already have a pos for this pos
      for(uint32_t j=0; j < root->amp_positions.size(); j++){
        position here_pos = root->amp_positions[j];
        if(here_pos.pos == tmp_pos[i].pos){
          found = true;
          root->amp_positions[j].depth += tmp_pos[i].depth;
          std::vector<allele> new_alleles = add_allele_vectors(tmp_pos[i].alleles, root->amp_positions[j].alleles);
          root->amp_positions[j].alleles = new_alleles;
          break;
        }
      }
      //if we've never seen this pos for this haplotype, push a new one
      if (!found){
        position tmp;
        tmp.depth = add_pos.depth;
        tmp.alleles = add_pos.alleles;
        tmp.pos = add_pos.pos;
        root->amp_positions.push_back(tmp);
      }
    }
   return; 
  }
}

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
bool doEnvelop(Interval i1, Interval i2) {
  if (i1.low <= i2.low && i1.high >= i2.high) return true;
  return false;
}

// The main function that searches an interval i in a given
// Interval Tree.
bool IntervalTree::envelopSearch(ITNode *root, Interval i) {
  // Base Case, tree is empty
  // std::cout << root->data->low << ":" << root->data->high << std::endl;
  if (root == NULL) return false;

  // If given interval overlaps with root
  if (doEnvelop(*(root->data), i)) return true;

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // be enveloped by an amplicon in left subtree
  if (root->left != NULL && root->left->max >= i.high)
    return envelopSearch(root->left, i);

  // Else interval can only be enveloped by amplicon in right subtree
  return envelopSearch(root->right, i);
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

/*
// Simple access functions to retrieve node's interval data
Interval ITNode::getData()const{
return data;
}
// Simple access functions to retrieve node's left child
ITNode ITNode::getLeft()const{
return left;
}
// Simple access functions to retrieve node's right child
ITNode ITNode::getRight()const{
return right;
}
// Simple access functions to set node's left child
void ITNode::setLeft(ITNode *node){
left = node;
}
// Simple access functions to set node's right child
void ITNode::setRight(ITNode *node){
right = node;
}

int main()
{
Interval ints[6] = {Interval(15, 20), Interval(30, 10), Interval(17, 19),
Interval(5, 20), Interval(12, 15), Interval(30, 40)}; int n = sizeof(ints) /
sizeof(ints[0]); IntervalTree tree = IntervalTree(); cout << "Hello World" <<
endl;
// populate interval tree
for (int i = 0; i < n; i++)
{
tree.insert(ints[i]);
}

tree.inOrder();
Interval queries[4] = {Interval(15, 20), Interval(9, 30), Interval(31, 38),
Interval(7, 22)}; int num_tests = sizeof(queries) / sizeof(queries[0]); for (int
i = 0; i < num_tests; i++)
{
cout << "Does " << queries[i].low << ":" << queries[i].high << " Overlap? " <<
tree.overlapSearch(queries[i]) << endl;
}
return 0;
}
*/
