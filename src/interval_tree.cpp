#include "interval_tree.h"
#include <chrono> 
using namespace std::chrono;
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


void IntervalTree::combine_haplotypes(ITNode *root){
  if (root==NULL) return;
  for(uint32_t i=0; i < root->amp_positions.size(); i++){
    int exists = check_position_exists(root->amp_positions[i].pos, variants);
    //does exist
    variants[root->amp_positions[i].pos].depth += root->amp_positions[i].depth;
    std::vector<allele> new_alleles = add_allele_vectors(root->amp_positions[i].alleles, variants[exists].alleles);
    variants[root->amp_positions[i].pos].alleles = new_alleles;
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

void IntervalTree::add_read_variants(uint32_t *cigar, uint32_t start_pos, uint32_t nlength, uint8_t *seq, uint8_t *aux, uint8_t* qual) {
  uint32_t consumed_query = 0;
  uint32_t consumed_ref = 0;
  std::vector<uint32_t> sc_positions; //sc positions            
  std::vector<uint32_t> ignore_sequence; //positions in query that are insertions                                  
  //first we handle insertion from the cigar
  std::vector<uint32_t> useful;
  std::vector<uint8_t> sequence;
  std::vector<uint8_t> quality;

  uint8_t nt = 0;
  uint32_t length = 0;
  uint32_t ll = 0;
  uint32_t op;
  for(uint32_t i=0; i < nlength; i++){
    op = bam_cigar_op(cigar[i]);
    if (bam_cigar_type(op) & 1){
      length = bam_cigar_oplen(cigar[i]);
      for(uint32_t k=0; k < length; k++){
        useful.push_back(ll+k);
      }
      ll += length;
    }
  }
  for(uint32_t k=0; k < useful.size(); k++){
    nt = seq_nt16_str[bam_seqi(seq, useful[k])];
    sequence.push_back(nt);
    quality.push_back(qual[useful[k]]);;
  }
  useful.clear();

  for(uint32_t j=0; j < nlength; j++){
    op = bam_cigar_op(cigar[j]);
    uint32_t oplen = bam_cigar_oplen(cigar[j]);
    if(op == 1){
      for(uint32_t k=0; k < oplen; k++){
        std::ostringstream convert;
        //convert data type to get the characters
        ignore_sequence.push_back(k+consumed_ref+start_pos);
        convert << sequence[k+consumed_query]; //this is throwing error
        std::string nuc = "+" + convert.str();
        //check if this position exists
        int  exists = check_position_exists(start_pos+consumed_ref, variants);
        if (exists != -1) {
          variants[exists].update_alleles(nuc,1, quality[k+consumed_query]);  
        } else {
          //add position to vector
          position add_pos;
          add_pos.pos = start_pos+consumed_ref; //don't add a position
          add_pos.update_alleles(nuc, 1, quality[k+consumed_query]);            
          variants.push_back(add_pos);
        }
      }
      consumed_query += oplen;
      continue;        
    }
    //if we don't consume both query and reference
    if (!(bam_cigar_type(op) & 1) || !(bam_cigar_type(op) & 2)){
      if (op != 2){
        for(uint32_t k=0; k < oplen; k++) {
          //convert data type to get the characters
          sc_positions.push_back(k+consumed_query);
        }
      }
    }
    //consumes query
    if (bam_cigar_type(op) & 1){
      consumed_query += oplen;
    } 
    //consumes ref
    if (bam_cigar_type(op) & 2){
      consumed_ref += oplen;
    } 
  }
  //we will use the aux tag to handle deletions
  bool deletion = false;
  uint32_t current_pos = start_pos;
  std::vector<uint32_t> deletion_positions; //the aux tag does NOT recognize insertions
  std::string gather_digits;     
  std::string deleted_char;    
  uint32_t last_char = 0;
  uint32_t j = 0;
  do {
    char character = (char)aux[j];
    if (character == '^'){
      current_pos += std::stoi(gather_digits);
      gather_digits = "";
      deletion = true;
    } else if (isdigit(character) && deletion) {
      deleted_char = "";
      deletion = false;
    } else if (isalpha(character) && deletion) {
      int exists = check_position_exists(current_pos, variants);
      if (exists != -1) {
        variants[exists].update_alleles("-", 1, 0);  
      } else {
        //add position to vector
        position add_pos;
        add_pos.pos = current_pos; //don't add a position
        add_pos.update_alleles("-", 1, 0);            
        variants.push_back(add_pos);
      } 
      deletion_positions.push_back(current_pos);
      current_pos += 1;
      deleted_char += character;
      deletion = true;    
    } else if (isdigit(character) && !deletion) {
      if(last_char > 0){
        current_pos += last_char;
        last_char = 0;
      } 
      gather_digits += character;
    } else if (isalpha(character) && !deletion) {
      last_char += 1;
      if(gather_digits.size() > 0){
        current_pos += std::stoi(gather_digits);
      }
      gather_digits = "";
    }
    j++; 
  } while(aux[j] != '\0');
  //now that we know where the insertions and deletions are, let's just iterate the query sequence and add it in, skipping problem positions
  current_pos = start_pos;
  //j is relative to the sequence and current pos to the reference
  //auto stop_again = high_resolution_clock::now();
  for(uint32_t j=0; j < sequence.size(); j++){
    std::vector<uint32_t>::iterator it = find(deletion_positions.begin(), deletion_positions.end(), current_pos);  
    if (it != deletion_positions.end()) {
      current_pos += 1;
      j -= 1;
      continue;
    }
    it = find(ignore_sequence.begin(), ignore_sequence.end(), current_pos);
    if (it != ignore_sequence.end()) {
      j += 1;
    }
    it = find(sc_positions.begin(), sc_positions.end(), j);
    if (it != sc_positions.end()){
      continue;
    }
    current_pos += 1;
    std::ostringstream convert;
    convert << sequence[j];
    std::string nuc = convert.str(); 
    variants[current_pos].update_alleles(nuc, 1, quality[j]);  
  }
  //auto stop_again_again = high_resolution_clock::now();
  //auto duration_again_again = duration_cast<microseconds>(stop_again_again - stop_again);
  //std::cerr << duration_again_again.count() << std::endl;
}

void IntervalTree::populate_variants(){
  for(uint32_t i=0; i <= max_pos; i++){
    position tmp;
    tmp.pos = i;
    variants.push_back(tmp);
  }
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
