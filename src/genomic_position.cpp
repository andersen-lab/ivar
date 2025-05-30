#include "genomic_position.h"
#include "interval_tree.h"
#include "saga.h"
#include <unordered_set>
#include <string>
#include <vector>

void get_amplicon_numbers(std::vector<amplicon_info> amplicons, std::vector<uint32_t> &amp_numbers){
  for(auto amp : amplicons){
    amp_numbers.push_back((uint32_t)amp.node->data->low);
  }
}

void set_amplicon_flag(std::vector<ITNode*> flagged_amplicons, std::vector<genomic_position> &global_positions){
  for(uint32_t i =0; i < global_positions.size(); i++){
    for(auto amp : global_positions[i].amplicons){
      ITNode* tmp = amp.node;
      bool exists = std::find(flagged_amplicons.begin(), flagged_amplicons.end(), tmp) != flagged_amplicons.end();
      if(exists){
        global_positions[i].amp_flux = true;
      }
    }
  }
}

void genomic_position::update_alleles(std::string nt, uint32_t qual){
  //check if in allele vector
  int exists = check_allele_exists(nt, alleles);
  //allele does not exist
  if (exists == -1){
    allele tmp;
    tmp.mean_qual = qual;
    tmp.depth = 1;
    tmp.nuc = nt;
    alleles.push_back(tmp);
    return;
  }
  alleles[exists].mean_qual += qual;
  alleles[exists].depth += 1;
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
    bool is_del = final_bases[i].find('-') != std::string::npos;
    bool is_ins = final_bases[i].find('+') != std::string::npos;

    uint32_t size = global_positions.size();
    while(size < final_positions[i]+1){
      genomic_position tmp;
      tmp.gapped_depth = 0;
      tmp.depth = 0;
      tmp.pos = global_positions.size();
      global_positions.push_back(tmp);
      size = global_positions.size();
    }
    uint32_t pos = final_positions[i];
    if(global_positions.size() >= pos){
      global_positions[final_positions[i]].update_alleles(final_bases[i], final_qualities[i]);
      if(is_del && !is_ins) global_positions[final_positions[i]].gapped_depth += 1;
      if(!is_del && !is_ins){
        global_positions[final_positions[i]].depth += 1;
        global_positions[final_positions[i]].gapped_depth += 1;
      }
    }
  }
}

void assign_read(ITNode *node, std::vector<uint32_t> final_positions, std::vector<std::string> final_bases, std::vector<uint32_t> final_qualities, std::vector<genomic_position> &global_positions) {
  for(uint32_t i=0; i < final_positions.size(); i++){
    uint32_t pos = final_positions[i];
    //for this position, iterate the amplicons associated
    bool found = false;
    bool is_del = final_bases[i].find('-') != std::string::npos;
    bool is_ins = final_bases[i].find('+') != std::string::npos;

    for(uint32_t j=0; j < global_positions[pos].amplicons.size(); j++){
      amplicon_info &amp = global_positions[pos].amplicons[j];
      found = node_compare(node, amp.node);
      if(found){
        amp.update_alleles(final_bases[i], final_qualities[i]);
        if (!is_del && !is_ins) {
          amp.amp_depth += 1;
          global_positions[pos].gapped_depth += 1;
          global_positions[pos].depth += 1;
        }
        if(!is_ins && is_del){
          amp.amp_depth_gapped += 1;
          global_positions[pos].gapped_depth += 1;
        }
        break;
      }
    }
    if(!found){
      //declare a new associated amplicon
      amplicon_info amp;
      amp.node = node;
      amp.amp_depth = 0;
      amp.amp_depth_gapped = 0;
      if (!is_del && !is_ins) {
        amp.amp_depth += 1;
        global_positions[pos].gapped_depth += 1;
        global_positions[pos].depth += 1;
      }
      if(!is_ins && is_del){
        amp.amp_depth_gapped += 1;
        global_positions[pos].gapped_depth += 1;
      }
      amp.amp_alleles = populate_basic_alleles();
      amp.update_alleles(final_bases[i], final_qualities[i]);
      global_positions[pos].amplicons.push_back(amp);
    }
  }
}

void collect_allele_frequencies(std::vector<amplicon_info> amplicons, std::unordered_map<std::string, std::vector<double>> &allele_frequencies) {
  for (auto amp : amplicons) {
    for (auto al : amp.amp_alleles) {
      if(al.depth == 0) continue;
      double freq = (double)al.depth / (double)amp.amp_depth;
      allele_frequencies[al.nuc].push_back(freq);
    }
  }
}

void collect_allele_stats(const std::vector<amplicon_info> &amplicons, std::unordered_map<std::string, std::vector<double>> &allele_frequencies, std::unordered_map<std::string, std::vector<uint32_t>> &depth_map, uint32_t min_depth, uint8_t min_qual){
  for (const auto &amp : amplicons) {
    for (const auto &al : amp.amp_alleles) {
      if (al.depth < min_depth || al.mean_qual < min_qual) continue;
      allele_frequencies[al.nuc].push_back(static_cast<double>(al.depth) / amp.amp_depth);
      depth_map[al.nuc].push_back(al.depth);
    }
  }
}
std::vector<ITNode*>  calculate_amplicon_variation(std::vector<genomic_position> &global_positions, uint32_t min_depth, uint8_t min_qual){
  std::vector<ITNode*> flagged_amplicons;
  std::unordered_map<std::string, std::vector<double>> allele_frequencies;
  std::unordered_map<std::string, std::vector<uint32_t>> allele_depths;
  std::unordered_set<ITNode*> seen_amplicons;
  for(uint32_t i=0; i < global_positions.size(); i++){
    if(global_positions[i].amplicons.size() > 0){
      allele_frequencies.clear();
      allele_depths.clear();
      collect_allele_stats(global_positions[i].amplicons, allele_frequencies, allele_depths, min_depth, min_qual);
      for (const auto &[key, values] : allele_frequencies) {
        double std = calculate_standard_deviation_weighted(values, allele_depths[key]);
        if(std > 0.03){
          global_positions[i].flux = true;
          //add the standard dev to the allele value
          for(auto &a : global_positions[i].alleles){
            if(a.nuc == key) a.stddev= std;
          }

          //add all amps to the flagged amps vec
          for(auto amp : global_positions[i].amplicons){
            ITNode* tmp = amp.node;
            if (seen_amplicons.insert(tmp).second) {
              flagged_amplicons.push_back(tmp);
            }
          }
        }
      }
    }
  }
  return(flagged_amplicons);
}

void add_allele_vectors(std::vector<allele> &alleles, const std::vector<allele> &amp_alleles){
  //build a hash map for fast lookup of existing alleles by nucleotide
  std::unordered_map<std::string, allele*> allele_map;
  for (auto &al : alleles) {
    allele_map[al.nuc] = &al;
  }
 //merge amplicon alleles into the map (and original vector)
  for (const auto &amp_al : amp_alleles) {
    auto it = allele_map.find(amp_al.nuc);
    if (it != allele_map.end()) {
      //update existing allele
      it->second->depth += amp_al.depth;
      it->second->mean_qual += amp_al.mean_qual;
    } else {
      //insert new allele and update map
      alleles.push_back(amp_al);
      allele_map[amp_al.nuc] = &alleles.back();
    }
  }
}

void combine_haplotypes(std::vector<genomic_position> &global_positions) {
  for(auto &position : global_positions){
    for(auto &amp : position.amplicons){
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
    tmp.depth = 0;
    tmp.gapped_depth = 0;
    positions.push_back(tmp);
  }
}
