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
  if (root==NULL) return;
  for(uint32_t i=0; i < root->amp_positions.size(); i++){
    if(root->amp_positions[i].depth == 0){
      continue;
    }
    variants[root->amp_positions[i].pos].depth += root->amp_positions[i].depth;
    std::vector<allele> new_alleles = add_allele_vectors(root->amp_positions[i].alleles, variants[root->amp_positions[i].pos].alleles);
    variants[root->amp_positions[i].pos].alleles = new_alleles;
    variants[root->amp_positions[i].pos].amplicon_numbers.push_back(counter);
    for(uint32_t j=0; j < root->amp_positions[i].alleles.size(); j++){
      if (variants[root->amp_positions[i].pos].amplicon_frequencies.find(root->amp_positions[i].alleles[j].nuc) != variants[root->amp_positions[i].pos].amplicon_frequencies.end()) {
        variants[root->amp_positions[i].pos].amplicon_frequencies[root->amp_positions[i].alleles[j].nuc].push_back((double)root->amp_positions[i].alleles[j].depth / (double)root->amp_positions[i].depth);
      } else{
        variants[root->amp_positions[i].pos].amplicon_frequencies[root->amp_positions[i].alleles[j].nuc] = {(double)root->amp_positions[i].alleles[j].depth / (double)root->amp_positions[i].depth};
      }
    }
  }
  counter += 1;
  combine_haplotypes(root->right, counter);
}

void IntervalTree::assign_read_amplicon(ITNode *root, uint32_t amp_start, std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities, uint8_t min_qual){
  if (root==NULL) return;
  if((uint32_t)root->data->low == amp_start){
    for(uint32_t i=0; i < positions.size(); i++){
      for(uint32_t j=0; j < root->amp_positions.size(); j++){
        if(positions[i] == root->amp_positions[j].pos){
          if(qualities[i] >= (uint32_t) min_qual){
            root->amp_positions[j].update_alleles(bases[i], 1, qualities[i]);
          }
        }
      }
    }
  }
  assign_read_amplicon(root->right, amp_start, positions, bases, qualities, min_qual);
}

void IntervalTree::find_read_amplicon(ITNode *root, uint32_t lower, uint32_t upper, bool &found, std::string read_name, uint32_t &amp_start, uint32_t &amp_dist){
  //read name here is for TEST
  if (root==NULL) return;
  //if ((uint32_t)root->data->low > upper) return;
  if(((uint32_t)root->data->low <= lower) && (upper <= (uint32_t)root->data->high)){
    //describes how far the ends of this are from the start/end of the amplicon
    uint32_t dist = (lower - root->data->low) + (root->data->high - upper);
    if(dist < amp_dist) { 
      amp_dist = dist;
      amp_start = root->data->low;
    }
    found = true;
  }
  find_read_amplicon(root->right, lower, upper, found, read_name, amp_start, amp_dist);
}

void IntervalTree::amplicon_position_pop(ITNode *root){
  if (root==NULL) return;
  for(uint32_t i=root->data->low; i < (uint32_t)root->data->high; i++){
    position add_pos;
    add_pos.pos = i;
    add_pos.alleles = populate_basic_alleles();
    root->amp_positions.push_back(add_pos);
  }
  amplicon_position_pop(root->right); 
}

void IntervalTree::detect_position_amplicons(ITNode *root, uint32_t find_position, uint32_t &counter, std::vector<uint32_t> &overlaps){
  if (root==NULL) return;
  if(((uint32_t)root->data->low < find_position) && (find_position < (uint32_t)root->data->high)){
    std::cerr << root->data->low << " " << find_position << " " << root->data->high << std::endl;
    overlaps.push_back(counter);
  }
  counter += 1;
  detect_position_amplicons(root->right, find_position, counter, overlaps);
}

void IntervalTree::detect_primer_issues(ITNode *root, uint32_t find_position){
  if (root==NULL) return;
  if (find_position < (uint32_t)root->data->low) return;
  if(((uint32_t)root->data->low == find_position) || (find_position == (uint32_t)root->data->high)){
    std::vector<uint32_t> tmp;
    tmp.push_back((uint32_t)root->data->low);
    tmp.push_back((uint32_t)root->data->high);
    overlaps.push_back(tmp);
    return;
  }
  detect_primer_issues(root->right, find_position);
}


void IntervalTree::detect_abberations(ITNode *root, uint32_t find_position){
  if (root==NULL) return;
    for(uint32_t i=0; i < root->amp_positions.size(); i++){
      if(find_position == root->amp_positions[i].pos){
        test_flux.push_back(root->amp_positions[i]);
        test_test.push_back(root->data->low);
        break;
      }
    }
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
    variants.push_back(tmp);
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
