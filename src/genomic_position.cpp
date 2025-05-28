#include "genomic_position.h"
#include "interval_tree.h"
#include <string>
#include <vector>

void genomic_position::update_alleles(std::string nt, uint32_t qual){
  //update overall positions depth
  if(nt.find("+") == std::string::npos){
    depth += 1;
  }
  //check if in allele vector
  int exists = check_allele_exists(nt, alleles);
  //allele does not exist
  if (exists == -1){
    allele tmp;
    tmp.mean_qual = qual;
    tmp.depth = 1;
    tmp.nuc = nt;
    alleles.push_back(tmp);
  } else {
    alleles[exists].mean_qual += qual;
    alleles[exists].depth += 1;
  }
}

void amplicon_info::update_alleles(std::string allele, uint32_t qual){
  for(uint32_t i=0; i < amp_alleles.size(); i++){
    if(amp_alleles[i].nuc == allele){
      amp_alleles[i].depth += 1;
      amp_alleles[i].mean_qual += qual;
    }
  }
}

void add_variants(std::vector<uint32_t> final_positions, std::vector<std::string> final_bases, std::vector<uint32_t> final_qualities, std::vector<genomic_position> &global_positions) {
  for(uint32_t i=0; i < final_positions.size(); i++){
    uint32_t size = global_positions.size();
    while(size-1 < final_positions[i]){
      genomic_position tmp;
      tmp.pos = final_positions[i];
      global_positions.push_back(tmp);
      size = global_positions.size();
    }
    uint32_t pos = final_positions[i];
    if(global_positions.size() >= pos){
      global_positions[final_positions[i]].update_alleles(final_bases[i], final_qualities[i]);
    }
  }
}

void assign_read(ITNode *node, std::vector<uint32_t> final_positions, std::vector<std::string> final_bases, std::vector<uint32_t> final_qualities, std::vector<genomic_position> &global_positions) {
  for(uint32_t i=0; i < final_positions.size(); i++){
    uint32_t pos = final_positions[i];
    //for this position, iterate the amplicons associated
    bool found = false;
    for(uint32_t j=0; j < global_positions[pos].amplicons.size(); j++){
      amplicon_info &amp = global_positions[pos].amplicons[j];
      found = node_compare(node, amp.node);
      if(found){
        amp.update_alleles(final_bases[i], final_qualities[i]);
        amp.amp_depth += 1;
        break;
      }
    }

    if(!found){
      //declare a new associated amplicon
      amplicon_info amp;
      amp.node = node;
      amp.amp_depth = 1;
      amp.amp_alleles = populate_basic_alleles();
      amp.update_alleles(final_bases[i], final_qualities[i]);
      global_positions[pos].amplicons.push_back(amp);
    }
  }
}

void calculate_amplicon_variation(std::vector<genomic_position> &global_positions){

}

void add_allele_vectors(std::vector<allele> &alleles, std::vector<allele> amp_alleles){
  for(uint32_t i=0; i < amp_alleles.size(); i++){
    bool found = false;
    for(uint32_t j=0; j < alleles.size(); j++){
      if(amp_alleles[i].nuc == alleles[j].nuc){
        found = true;
        alleles[j].depth += amp_alleles[i].depth;
        alleles[j].mean_qual += amp_alleles[i].mean_qual;
        break;
      }
    }
    //this allele hasn't been seen in the global yet
    if(!found){
      allele tmp;
      tmp.nuc = amp_alleles[i].nuc;
      tmp.depth = amp_alleles[i].depth;
      tmp.mean_qual = amp_alleles[i].mean_qual;
    }
  }
}

void combine_haplotypes(std::vector<genomic_position> &global_positions) {
  for(auto &position : global_positions){
    for(auto &amp : position.amplicons){
      position.depth += amp.amp_depth;
      add_allele_vectors(position.alleles, amp.amp_alleles);
    }
  }
}

int check_position_exists(uint32_t p, std::vector<genomic_position> positions) {
  for (uint32_t i=0; i < positions.size(); i++) {
    if (p == positions[i].pos) {
      return((int)i);
    }
  }
  return(-1);
}

void populate_positions(std::vector<genomic_position> &positions, uint32_t max_position){
  for(uint32_t i=0; i < max_position; i++) {
    genomic_position tmp;
    tmp.pos = i;
    positions.push_back(tmp);
  }
}
