#include "genomic_position.h"
#include "interval_tree.h"
#include "saga.h"
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

void amplicon_info::update_alleles(std::string allele_val, uint32_t qual){
  bool found = false;
  for(uint32_t i=0; i < amp_alleles.size(); i++){
    if(amp_alleles[i].nuc == allele_val){
      amp_alleles[i].depth += 1;
      amp_alleles[i].mean_qual += qual;
      found = true;
      break;
    }
  }
  //add a new allele if this one doesn't exist
  if(!found){
    allele tmp;
    tmp.nuc = allele_val;
    tmp.mean_qual = qual;
    tmp.depth = 1;
    amp_alleles.push_back(tmp);
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
        size_t pos = final_bases[i].find('-');
        if (pos == std::string::npos) {
          amp.amp_depth += 1;
        }
        amp.amp_depth_gapped += 1;
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

std::unordered_map<std::string, std::vector<uint32_t>> collect_allele_depths(std::vector<amplicon_info> amplicons) {
  std::unordered_map<std::string, std::vector<uint32_t>> depth_map;
  for (auto amp : amplicons) {
    for (auto al : amp.amp_alleles) {
      if(al.depth == 0) continue;
      depth_map[al.nuc].push_back(al.depth);
    }
  }
  return depth_map;
}

std::unordered_map<std::string, std::vector<double>> collect_allele_frequencies(std::vector<amplicon_info> amplicons) {
  std::unordered_map<std::string, std::vector<double>> freq_map;
  for (auto amp : amplicons) {
    for (auto al : amp.amp_alleles) {
      if(al.depth == 0) continue;
      double freq = (double)al.depth / (double)amp.amp_depth;
      freq_map[al.nuc].push_back(freq);
    }
  }
  return freq_map;
}

std::vector<uint32_t>  calculate_amplicon_variation(std::vector<genomic_position> &global_positions){
  std::vector<uint32_t> flagged_positions;
  for(uint32_t i=0; i < global_positions.size(); i++){
    if(global_positions[i].amplicons.size() > 0){
      std::vector<amplicon_info> tmp = global_positions[i].amplicons;
      std::unordered_map<std::string, std::vector<double>> allele_frequencies = collect_allele_frequencies(tmp);
      std::unordered_map<std::string, std::vector<uint32_t>> allele_depths = collect_allele_depths(tmp);
      for (const auto &[key, values] : allele_frequencies) {
        double std = calculate_standard_deviation_weighted(values, allele_depths[key]);
        if(std > 0.03){
          global_positions[i].flux = true;
          bool exists = std::find(flagged_positions.begin(), flagged_positions.end(), global_positions[i].pos) != flagged_positions.end();
        }
      }
    }
  }
  return(flagged_positions);
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

uint32_t calculate_gapped_depth(std::vector<allele> alleles){
  uint32_t gapped_depth = 0;
  for(auto a : alleles){
    size_t pos = a.nuc.find('-');
    if (pos != std::string::npos) {
      gapped_depth += a.depth;
    }
  }
  return(gapped_depth);
}

void combine_haplotypes(std::vector<genomic_position> &global_positions) {
  for(auto &position : global_positions){
    for(auto &amp : position.amplicons){
      uint32_t gapped_depth = calculate_gapped_depth(amp.amp_alleles);
      position.depth += amp.amp_depth - gapped_depth;
      position.gapped_depth += amp.amp_depth;
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
