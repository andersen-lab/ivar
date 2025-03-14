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

void IntervalTree::combine_haplotypes(ITNode *root){
  if (root==NULL) return;
  for(uint32_t i=0; i < root->amp_positions.size(); i++){
    if(root->amp_positions[i].depth == 0){
      continue;
    }
    variants[root->amp_positions[i].pos].depth += root->amp_positions[i].depth;
    std::vector<allele> new_alleles = add_allele_vectors(root->amp_positions[i].alleles, variants[root->amp_positions[i].pos].alleles);
    variants[root->amp_positions[i].pos].alleles = new_alleles;
  }
  combine_haplotypes(root->right);
}

void IntervalTree::find_read_amplicon(ITNode *root, uint32_t lower, uint32_t upper, std::vector<uint32_t> positions, std::vector<std::string> bases, std::vector<uint32_t> qualities){
  if (root==NULL) return;
  //if ((uint32_t)root->data->low > upper) return;
  if(((uint32_t)root->data->low <= lower) && (upper <= (uint32_t)root->data->high)){
    for(uint32_t i=0; i < positions.size(); i++){
      for(uint32_t j=0; j < root->amp_positions.size(); j++){
        if(positions[i] == root->amp_positions[j].pos){
          root->amp_positions[j].update_alleles(bases[i], 1, qualities[i]);
        }
      }
    }
  }
  find_read_amplicon(root->right, lower, upper, positions, bases, qualities);
}

void IntervalTree::amplicon_position_pop(ITNode *root){
  if (root==NULL) return;
  for(uint32_t i=root->data->low; i < root->data->high; i++){
    position add_pos;
    add_pos.pos = i;
    add_pos.alleles = populate_basic_alleles();
    root->amp_positions.push_back(add_pos);
  }
  amplicon_position_pop(root->right); 
}

void IntervalTree::detect_amplicon_overlaps(ITNode *root, uint32_t find_position){
  if (root==NULL) return;
  if (find_position < (uint32_t)root->data->low) return;
  if(((uint32_t)root->data->low < find_position) && (find_position < (uint32_t)root->data->high)){
    std::vector<uint32_t> tmp;
    tmp.push_back((uint32_t)root->data->low);
    tmp.push_back((uint32_t)root->data->high);
    overlaps.push_back(tmp);
  }
  detect_amplicon_overlaps(root->right, find_position);
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
  //unsure why this line is wrogn, but it is
  //if (find_position < (uint32_t)root->data->low) return;
    for(uint32_t i=0; i < root->amp_positions.size(); i++){
      if(find_position == root->amp_positions[i].pos){
        test_flux.push_back(root->amp_positions[i]);
        test_test.push_back(root->data->low);
        break;
      }
    }
  detect_abberations(root->right, find_position);

}

void IntervalTree::add_read_variants(uint32_t *cigar, uint32_t start_pos, uint32_t nlength, uint8_t *seq, uint8_t *aux, uint8_t* qual, std::string qname) {
  uint32_t consumed_query = 0;
  uint32_t consumed_ref = 0;
  std::vector<uint32_t> sc_positions; //sc positions            
  std::vector<uint32_t> ignore_sequence; //positions in query that are insertions                                  
  //first we handle insertion from the cigar
  std::vector<uint32_t> useful;
  std::vector<uint8_t> sequence;
  std::vector<uint8_t> quality;

  uint32_t quality_threshold=20;
  uint8_t nt = 0;
  uint32_t length = 0, ll=0, op=0;
  //auto start = high_resolution_clock::now();
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
  //TESTLINES
  uint32_t qual_threshold = 20;
  quality.reserve(useful.size());
  sequence.reserve(useful.size());
  for(uint32_t k=0; k < useful.size(); k++){
    if(qual[useful[k]] < qual_threshold){
      sequence.push_back('L');
      total_qual += qual[useful[k]];
      total_bases++;
    } else {
      nt = seq_nt16_str[bam_seqi(seq, useful[k])];
      sequence.push_back(nt);
    }
    quality.push_back(qual[useful[k]]);;
  }
  useful.clear();
  std::string nuc;
  for(uint32_t j=0; j < nlength; j++){
    op = bam_cigar_op(cigar[j]);
    uint32_t oplen = bam_cigar_oplen(cigar[j]);
    if(op == 1){
      std::ostringstream convert;
      uint32_t q = 0;
      for(uint32_t k=0; k < oplen; k++) {
        //convert data type to get the characters
        ignore_sequence.push_back(k+consumed_ref+start_pos);
        convert << sequence[k+consumed_query];
        q += quality[k+consumed_query];
      }
      nuc.clear();
      nuc = convert.str();

      int avg_q = (int)q/nuc.size();
      char ch = 'L';
      std::string nuc = "+" + convert.str();
      //check if this position exists
      int  exists = check_position_exists(start_pos+consumed_ref, variants);
      if (exists != -1 &&  nuc.find(ch) == std::string::npos) {
        variants[exists].update_alleles(nuc, 1, avg_q);  
      } else {
        //add position to vector
        position add_pos;
        add_pos.pos = start_pos+consumed_ref; //don't add a position
        if (nuc.find(ch) == std::string::npos) {
          add_pos.update_alleles(nuc, 1, avg_q);            
          variants.push_back(add_pos);
        }
      }
      consumed_query += oplen;
      continue;        
    }
    if (!(bam_cigar_type(op) & 1) && !(bam_cigar_type(op) & 2)){
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
  bool substitution = false;
  uint32_t current_pos = start_pos;
  std::vector<uint32_t> deletion_positions; //the aux tag does NOT recognize insertions
  std::vector<uint32_t> substitutions; //handle substitutions
  std::string gather_digits = "";
  std::string deleted_char;
  std::string sub_char;
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
      gather_digits += character;
    } else if (isalpha(character) && deletion) {
      variants[current_pos].update_alleles("-", 1, 0);  
      deletion_positions.push_back(current_pos);
      current_pos += 1;
      deleted_char += character;
      deletion = true;
    } else if (isdigit(character) && !deletion) {
      if(substitution){
        for(uint32_t z=1; z < sub_char.size(); z++){
          substitutions.push_back(current_pos + z);
        }
        substitution = false;
        sub_char.clear();
      }
      if(last_char > 0){
        current_pos += last_char;
        last_char = 0;
      }
      gather_digits += character;
    } else if (isalpha(character) && !deletion) {
      last_char += 1;
      substitution = true;
      sub_char += character;
      if(gather_digits.size() > 0){
        current_pos += std::stoi(gather_digits);
      }
      gather_digits = "";
    }
    j++;
  } while(aux[j] != '\0');
  //now that we know where the insertions and deletions are, let's just iterate the query sequence and add it in, skipping problem positions
  current_pos = start_pos;
  bool prev_insertion = false;
  std::ostringstream convert;
  nuc.clear();
  std::vector<uint32_t> seen_insertions;
  std::vector<uint32_t>::iterator i_it;
  std::string test = "";
  //j is relative to the sequence and current pos to the reference
  for(uint32_t j=0; j < sequence.size(); j++){
    if(qname == test){
      std::cerr << "j " << j << " cur pos " << current_pos << " " << sequence[j] << std::endl;
    }
    std::vector<uint32_t>::iterator it = find(deletion_positions.begin(), deletion_positions.end(), current_pos); 
    if (it != deletion_positions.end()) {
      current_pos += 1;
      j -= 1;
      continue;
    }
    it = find(ignore_sequence.begin(), ignore_sequence.end(), current_pos);
    i_it = find(seen_insertions.begin(), seen_insertions.end(), current_pos);
    //handle insertions
    if (it != ignore_sequence.end() && i_it == seen_insertions.end()) {
      std::ostringstream convert;
      convert << sequence[j];
      nuc += convert.str();
      seen_insertions.push_back(current_pos);
      current_pos += 1;
      prev_insertion = true;
      continue;
    } else if (it == ignore_sequence.end()){
      //insertion is over 
      if(prev_insertion) {
        current_pos -= nuc.size();
      }
      prev_insertion = false;
      nuc = "";
    }
    it = find(sc_positions.begin(), sc_positions.end(), j);
    if (it != sc_positions.end()){
      continue;
    }
    current_pos += 1;
    std::ostringstream convert;
    convert << sequence[j];
    std::string nuc = convert.str(); 
    if((uint32_t)quality[j] < quality_threshold){
      continue;
    }
    char ch = 'L';
    if (nuc.find(ch) == std::string::npos) {
      variants[current_pos].update_alleles(nuc, 1, quality[j]);  
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
