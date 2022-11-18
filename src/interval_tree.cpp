#include "interval_tree.h"
#include "clustering.h"
// Constructor for initializing an Interval Tree
IntervalTree::IntervalTree(){
  _root = NULL;
}

//clear all values associated with amplicons
void IntervalTree::clear(ITNode *root){
  /*
   * Loop through the interval tree and clear data from each node.
   */
  if(root == NULL){
    return;
  }
  root->read_count = 0;
  root->frequency.clear();
  root->haplotypes.clear();
  root->positions.clear();
  root->ranges.clear();
  root->final_positions.clear();
  root->final_haplotypes.clear(); 
  clear(root->right);
}

// A utility function to insert a new Interval Search Tree Node
// This is similar to BST Insert.  Here the low value of interval
// is used tomaintain BST property
void IntervalTree::insert(ITNode *root, Interval data){
  // Base case: Tree is empty, new node becomes root
  if(root == NULL){
    root = new ITNode(data);
    _root = root;
  } else {
    // Get low value of interval at root
    uint32_t l = root->data->low;
    // If root's low value is greater, then new interval goes to
    // left subtree
    if (data.low < l){
      if(!root->left){
	ITNode *tmpNode = new ITNode(data);
	//std::cout << data.low << ":" << data.high << "->insertLeft" << std::endl;
	root->left = tmpNode;
      } else {
	insert(root->left, data);
      }
    }
    else {
      if(!root->right){
	ITNode *tmpNode = new ITNode(data);
	//std::cout << data.low << ":" << data.high << "->insertRight" << std::endl;
	root->right = tmpNode;
      } else {
	insert(root->right, data);
      }
    }
  }
  // update max value of ancestor node
  if(root->max < data.high)
    root->max = data.high;
}


// A utility function to check if the 1st interval envelops the second
bool doEnvelop(Interval i1, Interval i2){
  if(i1.low <= i2.low && i1.high >= i2.high)
    return true;
  return false;
}


// The main function that searches an interval i in a given
// Interval Tree.
bool IntervalTree::envelopSearch(ITNode *root, Interval i){
  // Base Case, tree is empty
  //std::cout << root->data->low << ":" << root->data->high << std::endl;
  if (root == NULL) return false;

  // If given interval overlaps with root
  if (doEnvelop(*(root->data), i))
    return true;

  // If left child of root is present and max of left child is
  // greater than or equal to given interval, then i may
  // be enveloped by an amplicon in left subtree
  if (root->left != NULL && root->left->max >= i.high)
    return envelopSearch(root->left, i);

  // Else interval can only be enveloped by amplicon in right subtree
  return envelopSearch(root->right, i);
}

// A helper function for inorder traversal of the tree
void IntervalTree::inOrder(ITNode *root){
  if (root == NULL) return;
  inOrder(root->left);
  cout << "[" << root->data->low << ", " << root->data->high << "]"
       << " max = " << root->max << endl;
  inOrder(root->right);
}

bool primer_binding_mutation(std::vector<uint32_t> positions, uint32_t primer_start, uint32_t primer_end){
  /*
   * @param positions : vector containing all positions mutation occurs at
   * @param primer_start : the start pos of the primer  
   * @params primer_end : the end pos of the primer
   *
   * Function looks for a mutation in the primer binding region for side of the read
   * which hasn't been soft-clipped.
   */
  bool primer_mutation = false;
  for(uint32_t i = 0; i < positions.size(); i++){
    if((positions[i] < primer_end) && (positions[i] > primer_start)){
      primer_mutation = true;
      return(primer_mutation);  
    }
  }
  return(primer_mutation);
}

std::vector<std::pair<std::uint32_t, int>>  _trim_read_positions(std::vector<int> haplotypes, std::vector<uint32_t> positions, uint32_t lower_bound, uint32_t upper_bound){
  
  /* @param haplotypes : vector containing haplotype information
   * @param positions : vector containing the positions 
   * @param lower_bound : the lower bound of the amplicon
   * @param upper_bound : the upper bound of the amplicon
   * @param zipped : the zipped haplotypes and positions
   *
   * Function takes in the haplotypes and positions of an amplicon and removes any values that don't
   * fall within the amplicon. 
   */
  
  std::vector<int> new_haplotypes;
  std::vector<uint32_t> new_positions;
  for(uint32_t i=0; i < haplotypes.size(); i++){
    //in range of the amplicon
    if((positions[i] >= lower_bound) && (positions[i] <= upper_bound)){
      new_haplotypes.push_back(haplotypes[i]);
      new_positions.push_back(positions[i]);
    }
  }

  //in order to pass two values back we 'zip' them together
  std::vector<std::pair<std::uint32_t, int>> zipped;
  zip(new_haplotypes, new_positions, zipped);
  return(zipped);
}

void remove_low_quality_nts(ITNode *node, std::vector<position> all_positions){
  /*
   * Function takes in a single node and uses average quality information to eliminate positions
   * within haplotypes that are possible sequencing errors.
   */
  if (node == NULL) return;
  std::vector<allele> allele_storage; 
  std::string nt;
  int one = -1;
  for(uint32_t i = 0; i < node->positions.size(); i++){
    for(uint32_t x = 0; x < node->positions[i].size(); x++){
      //we don't care about pos alredy SC or match ref
      if(node->haplotypes[i][x] < 0){
        continue;
      }
      nt = decoded_nucs(node->haplotypes[i][x]);
      allele_storage = all_positions[node->positions[i][x]].ad;
      for (allele a : allele_storage){
        //find the proper nuc
        if(nt == a.nuc){
          //check the quality
          if(a.mean_qual+0 < 20){
            //if quality is low we pretend it's soft clipped
            node->haplotypes[i][x] = one;
          }
        }
      }
    }
  }
}

//TODO estimate the location of the read in relation to the total genome for performance gains
/*void(){

}*/

//traverse the tree and find the amplicon the read belongs within
void IntervalTree::find_amplicon_per_read(ITNode *root, uint32_t start, uint32_t end, 
    std::vector<int> haplotypes, std::vector<uint32_t> positions, bool reverse,
    std::vector<uint32_t> range, std::vector<position> &all_positions, bool primer_mutation, bam1_t *r){
  /*
   * @param root : IT node obj
   * @param start : left-most position of read
   * @param end : right-most position of read
   * @param haplotypes : the haplotypes (NT) to be placed on an amplicon
   * @param positions : the positons that go with the haplotypes 
   * @param reverse : whether or not a read is reversed
   * @param range : marks the start and end of the read
   * @param all_positions : vector containing stored info for all positions (ie. alleles, depths etc.)
   * @param primer_mutation : is there a mutation in the primer binding region
   *
   * Takes in an interval tree, and a read start and end pos and adds 
   * the haplotype, positions, and frequencies to to internval tree object.
   */
   
  if (root == NULL) return;
  //yes karthik, I know this is poorly written
  if(reverse){
    if((root->data->high_inner+1 - 10 < end) && (root->data->high_inner + 10 + 1 > end)){
      root->reverse += 1;
      if(primer_mutation){
        root->mut_reverse += 1;
      }
      //we go through this additional step where we chop off positions not within the amplicon
      std::vector<std::pair<std::uint32_t, int>> zipped = _trim_read_positions(haplotypes, positions, 
          root->data->low, root->data->high);
        if(zipped.size() == 0){return;}
        //unzip the newly modifed haplotypes
        unzip(zipped, haplotypes, positions);
        if(zipped.size() == 0){return;}
        root->haplotypes.push_back(haplotypes);
        root->positions.push_back(positions);
        root->ranges.push_back(range);
        root->read_count += 1;
       return;
    }
  }else{
    if((root->data->low_inner - 10 - 1 < start) && (root->data->low_inner + 10 - 1 >  start)){
      root->forward += 1;
      if(primer_mutation){
        root->mut_forward += 1;
      }     
     std::vector<std::pair<std::uint32_t, int>> zipped = _trim_read_positions(haplotypes, positions, 
        root->data->low, root->data->high);
      if(zipped.size() == 0){return;}
      //unzip the newly modifed haplotypes
      unzip(zipped, haplotypes, positions);
      root->haplotypes.push_back(haplotypes);
      root->positions.push_back(positions);
      root->ranges.push_back(range);
      //always increment the read count, even if the read mathes the ref perfectly
      root->read_count += 1;
      return;
   }
  }
  //make sure we haven't exceeded the possible range for amplicons for this read
  find_amplicon_per_read(root->right, start, end, haplotypes, positions, reverse, range, all_positions, primer_mutation, r);
  
}

//use this to iterate through the tree, returning each node
ITNode* IntervalTree::iterate_nodes(ITNode *root){
  /*
   * @param root : node of the tree
   *
   * Function takes an IT node and returns it, much like a get function.
   */
  if(root == NULL) return(root);
  return(root);
}

//use this to dump amplicon summary data to a tab seperated text file
void IntervalTree::dump_amplicon_summary(ITNode *root, std::string filename){
  /*
   * @param root : node of the interval tree
   * @param filename : full path to file where amplicon information is stored
   *
   * Opens tab seperated text file and writes information related to amplicon level
   * haplotypes and frequencies.
   */
  ofstream file;
  if (root == NULL) return;
  if ((root->read_count != 0) && (root->final_haplotypes.size() > 0)){
    //open the amplicon info file
    file.open(filename, ios_base::app);
    //set the file headers
    //file << "lower_primer\tupper_primer\tread_count\tpositions\tfrequencies\thaplotypes\n";
    file << root->data->low << "\t";
    file << root->data->high << "\t";
    file << root->read_count << "\t";  
    //std::cout << root->data->low << std::endl;
    int count = 0;
    std::string positions;
    for(uint32_t x:root->final_positions){
      if(count !=0){ positions += "_";}
      //std::cout << x << " ";
      count += 1;
      positions += std::to_string(x);
    }
    //std::cout << "\n";
    file << positions << "\t"; 
    std::string frequency;
    count = 0;

    for(float i:root->frequency){
      if(count !=0){frequency += "_";}
      count += 1;
      //std::cout << i << " ";
      frequency += std::to_string(i);
    }
    //std::cout << "\n";
    file << frequency << "\t";
    
    std::string haplotypes;
    std::string tmp;
    count = 0;
    int count_tmp = 0;
    //something weird going on with this
    for(std::vector<int> hap:root->final_haplotypes){
      tmp.clear();
      count_tmp = 0;
      for(int h:hap){
        //std::cout << h << " ";
        if(count_tmp != 0){tmp += "_";}
        tmp += std::to_string(h);
        count_tmp += 1;
      }
      //std::cout << "\n";
      if(count != 0){haplotypes += ";";}
      haplotypes += tmp;
      count += 1;
    }
    file << haplotypes << "\t";
    file << root->final_haplotypes.size() << std::endl;
    file.close();
  }
  dump_amplicon_summary(root->right, filename);
  
}

//use this to print out summayr of unique haplotypes and frequencies
void IntervalTree::print_amplicon_summary(ITNode *root){
  if (root == NULL) return;
  if ((root->read_count != 0) && (root->final_positions.size() > 0)){
    std::cout << "\nLow:" << root->data->low << " High: " << root->data->high << " ReadCount: " << root->read_count << std::endl;
    //print position info
    for(uint32_t x:root->final_positions){
        std::cout << x << ", ";
    }
    std::cout << "\n";  
    for(float i:root->frequency){
      std::cout << i << ", ";
    }
    std::cout << "\n";
    //print haplotype info
    for(std::vector<int> hap:root->final_haplotypes){
      for(int h:hap){
        std::cout << h << ", ";
      }
      std::cout << "\n";
    }
  }
  print_amplicon_summary(root->right);
}

// A stand-alone function to create a tree containing the coordinates of each amplicon
// based on user-specified primer pairs
IntervalTree populate_amplicons(std::string pair_info_file, std::vector<primer> &primers){
  int amplicon_start = -1;
  int amplicon_end = -1;
  int amplicon_start_inner = -1;
  int amplicon_end_inner = -1;

  IntervalTree tree = IntervalTree();
  populate_pair_indices(primers, pair_info_file);
  for (auto & p : primers) {
    if (p.get_strand() == '+'){
      if (p.get_pair_indice() != -1){
	      amplicon_start = p.get_start();
	      amplicon_end = primers[p.get_pair_indice()].get_end() + 1;
        amplicon_start_inner = p.get_end() +1;
        amplicon_end_inner = primers[p.get_pair_indice()].get_start();
	      tree.insert(Interval(amplicon_start, amplicon_end, amplicon_start_inner, amplicon_end_inner));
	    }
    }
  }
  return tree;
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
Interval ints[6] = {Interval(15, 20), Interval(30, 10), Interval(17, 19), Interval(5, 20), Interval(12, 15), Interval(30, 40)};
int n = sizeof(ints) / sizeof(ints[0]);
IntervalTree tree = IntervalTree();
cout << "Hello World" << endl;
// populate interval tree
for (int i = 0; i < n; i++)
{
tree.insert(ints[i]);
}

tree.inOrder();
Interval queries[4] = {Interval(15, 20), Interval(9, 30), Interval(31, 38), Interval(7, 22)};
int num_tests = sizeof(queries) / sizeof(queries[0]);
for (int i = 0; i < num_tests; i++)
{
cout << "Does " << queries[i].low << ":" << queries[i].high << " Overlap? " << tree.overlapSearch(queries[i]) << endl;
}
return 0;
}
*/
