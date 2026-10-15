#include "call_consensus_clustering.h"
#include "gmm.h"
#include "saga.h"
#include "allele_functions.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>

consensus_sequence build_majority_consensus(const std::vector<variant> &variants, double threshold){
  uint32_t max_position = 0;
  for(const auto &x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  }
  consensus_sequence cs(max_position);
  //a majority consensus votes over every record, so don't filter on consensus_numbers
  //the way assign_variants_position does
  for(uint32_t i=0; i < variants.size(); i++){
    uint32_t span = 1;
    if(variants[i].nuc.find('-') != std::string::npos){
      std::string nuc = variants[i].nuc;
      nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
      span = (uint32_t)nuc.size();
    }
    for(uint32_t z=0; z < span && variants[i].position + z <= max_position; z++){
      cs.add_variant(variants[i].position + z, variants[i]);
    }
  }
  cs.process_variant_assignments();
  cs.get_majority_consensus(threshold);
  return cs;
}

void call_majority_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold){
  std::cerr << "in majority consensus call" << std::endl;
  consensus_sequence cs = build_majority_consensus(variants, default_threshold);
  cs.set_seq_name(clustering_file + "_" + std::to_string(default_threshold) + "_threshold");
  std::string consensus_filename = clustering_file + "_threshold.fa";
  std::ofstream(consensus_filename, std::ios::trunc); //start each run from an empty file, since write_consensus_to_file appends
  cs.write_consensus_to_file(consensus_filename);
}

void assign_variants_position(std::vector<variant> &variants, std::vector<consensus_sequence> &all_consensus_seqs){
  for(uint32_t i=0; i < variants.size(); i++){
    //deletions span multiple positions; register the record at every position
    //it covers, not just its own start, so downstream per-position lookups
    //(and the deletion-span walk in process_variant_assignments) can find it.
    uint32_t span = 1;
    if(variants[i].nuc.find('-') != std::string::npos){
      std::string nuc = variants[i].nuc;
      nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
      span = (uint32_t)nuc.size();
    }
    //no population credibly carries this allele, so the position is unresolved in every
    //genome, not just the ones it was assigned to
    if(variants[i].freq_outlier){
      for(uint32_t k=0; k < all_consensus_seqs.size(); k++){
        for(uint32_t z=0; z < span; z++){
          all_consensus_seqs[k].mark_ambiguous(variants[i].position + z);
        }
      }
      continue;
    }
    for(uint32_t j=0; j < variants[i].consensus_numbers.size(); j++){
      uint32_t k = variants[i].consensus_numbers[j];
      for(uint32_t z=0; z < span; z++){
        all_consensus_seqs[k].add_variant(variants[i].position + z, variants[i]);
      }
    }
    for(uint32_t j=0; j < variants[i].ambiguous_numbers.size(); j++){
      uint32_t k = variants[i].ambiguous_numbers[j];
      for(uint32_t z=0; z < span; z++){
        all_consensus_seqs[k].mark_ambiguous(variants[i].position + z);
      }
    }
  }
}

void consensus_sequence::get_consensus(){
  for(uint32_t i=0; i < variant_records.size(); i++){
    if(position_ambiguous[i]){
      continue; //multiple explanations for this peak - leave as N
    }
    //frequency fluctuates across amplicons - leave as N
    if(std::any_of(variant_records[i].begin(), variant_records[i].end(),
                   [](const variant &v){ return v.position_masked && !v.half_normal_upper; })){
      continue;
    }
    if(variant_records[i].size() == 0){
      continue;
    }
    std::vector<string> nucs; //all the possible nucs here?
    bool deletion_added = false;
    std::string insertion;
    double insertion_freq = -1.0; //-1 not 0: gapped_freq defaults to 0
    for(uint32_t j=0; j < variant_records[i].size(); j++){
      bool has_insertion = variant_records[i][j].nuc.find('+') != std::string::npos;
      if(has_insertion){
        const variant &v = variant_records[i][j];
        if(v.qual_flag || v.depth_flag || v.half_normal_lower) continue;
        std::string candidate = v.nuc;
        candidate.erase(std::remove(candidate.begin(), candidate.end(), '+'), candidate.end());
        //keep the best supported, ties broken on length then lexically
        if(v.gapped_freq > insertion_freq ||
           (v.gapped_freq == insertion_freq &&
            (candidate.size() > insertion.size() ||
             (candidate.size() == insertion.size() && candidate < insertion)))){
          insertion = candidate;
          insertion_freq = v.gapped_freq;
        }
        //an insertion rides along with the base called here, it is not an allele
        //of its own - letting "+..." reach nucs masks the position via gt2iupac
        continue;
      }

      //if we have one assigned to the upper half normal, we go with that and ignore the rest

      if(variant_records[i][j].half_normal_upper){
        //deletions are written as a gap, as in the normal path below
        if(variant_records[i][j].nuc.find('-') != std::string::npos){
          if(!deletion_added){
            nucs.push_back("-");
            deletion_added = true;
          }
        } else {
          nucs.push_back(variant_records[i][j].nuc);
        }
        continue;
      }
      if(variant_records[i][j].position_half_normal_upper){
        continue;
      }
      //if we have it assigned to the lower half normal, we ignore it
      if(variant_records[i][j].half_normal_lower){
        continue;
      }

      if(variant_records[i][j].qual_flag || variant_records[i][j].depth_flag){
        continue;
      }

      if(variant_records[i][j].assigned_deletion && !deletion_added){
        nucs.push_back("-");
        deletion_added = true;
      }

      if(variant_records[i][j].nuc.find('-') == std::string::npos){
        nucs.push_back(variant_records[i][j].nuc);
      }
    } 

    //we only have one variant at this position, so we can assign it to the consensus sequence
    if(nucs.size() == 1){
      sequence[i] = nucs[0] + insertion;
    } else if(nucs.size() > 1){
      //we have multiple variants and one is a deletion so we call an N
      bool has_deletion = std::any_of(nucs.begin(), nucs.end(), [](const std::string &s){ return s == "-"; });
      if(has_deletion){
        continue;
      }
      //we have multiple variants at this position, so we need to combine them into an IUPAC code
      char combined = nucs[0][0];
      for(uint32_t j=1; j < nucs.size(); j++){
        combined = gt2iupac(combined, nucs[j][0]);
      }
      sequence[i] = std::string(1, combined) + insertion;
    }
  }
}

void consensus_sequence::get_majority_consensus(double threshold){
  for(uint32_t i=0; i < variant_records.size(); i++){
    if(variant_records[i].size() == 0){
      continue;
    }

    std::vector<std::string> nucs;
    std::vector<double> freqs;
    std::string insertion;
    double insertion_freq = -1.0;

    for(uint32_t j=0; j < variant_records[i].size(); j++){
      if(variant_records[i][j].qual_flag || variant_records[i][j].depth_flag){
        continue;
      }

      bool has_insertion = variant_records[i][j].nuc.find('+') != std::string::npos;
      if(has_insertion){
        const variant &v = variant_records[i][j];
        //insertion frequency is not part of the substitution sum, so test it
        //against the threshold on its own rather than adding it to freqs[]
        if(v.gapped_freq >= threshold){
          std::string candidate = v.nuc;
          candidate.erase(std::remove(candidate.begin(), candidate.end(), '+'), candidate.end());
          //several insertion alleles can coexist here; keep the best supported
          //rather than the last one read, ties broken on length then lexically
          if(v.gapped_freq > insertion_freq ||
             (v.gapped_freq == insertion_freq &&
              (candidate.size() > insertion.size() ||
               (candidate.size() == insertion.size() && candidate < insertion)))){
            insertion = candidate;
            insertion_freq = v.gapped_freq;
          }
        }
        continue;
      }

      //read the deletion off the record itself, not the assigned_deletion stamp -
      //process_variant_assignments stamps every allele across a deletion's span, which
      //would let a sub-threshold deletion erase the real bases there
      if(variant_records[i][j].nuc.find('-') != std::string::npos){
        nucs.push_back("-");
      } else {
        nucs.push_back(variant_records[i][j].nuc);
      }
      freqs.push_back(variant_records[i][j].gapped_freq);
    }

    if(nucs.empty()) continue;

    //sort indices by frequency, descending
    std::vector<uint32_t> order(nucs.size());
    for(uint32_t j=0; j < order.size(); j++) order[j] = j;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b){ return freqs[a] > freqs[b]; });

    //walk alleles largest -> smallest, pulling in an entire tied-frequency
    //group at once, accumulating frequency until the running sum reaches or
    //exceeds threshold (or every allele has been included)
    std::vector<std::string> included;
    double cumulative = 0.0;
    uint32_t idx = 0;
    while(idx < order.size()){
      double group_freq = freqs[order[idx]];
      while(idx < order.size() && freqs[order[idx]] == group_freq){
        included.push_back(nucs[order[idx]]);
        cumulative += freqs[order[idx]];
        idx++;
      }
      if(cumulative >= threshold) break;
    }

    //every allele was included and the total still fell short of the threshold,
    //so no call here is supported - leave it as N
    if(cumulative < threshold) continue;

    if(included.size() == 1){
      sequence[i] = included[0] + insertion;
    } else {
      //multiple alleles needed to reach the threshold - call an ambiguity between them
      bool has_deletion = std::any_of(included.begin(), included.end(), [](const std::string &s){ return s == "-"; });
      if(has_deletion){
        continue; //deletion mixed with alleles can't be IUPAC-combined - leave as "N"
      }
      char combined = included[0][0];
      for(uint32_t j=1; j < included.size(); j++){
        combined = gt2iupac(combined, included[j][0]);
      }
      sequence[i] = std::string(1, combined) + insertion;
    }
  }
}

void consensus_sequence::process_variant_assignments(){
  //flag positions where deletion overlap, and where we have an upper half normal assignment
  for(uint32_t i=0; i < variant_records.size(); i++){
    //anything in the upper half normal wins the position, whatever order the records are in
    bool position_half_normal_upper = std::any_of(variant_records[i].begin(), variant_records[i].end(),
                                                  [](const variant &v){ return v.half_normal_upper; });
    bool is_del = false;
    uint32_t deletion_span = 0;
    for(const auto &v : variant_records[i]){
      if(v.nuc.find('-') == std::string::npos) continue;
      //copies sit at every covered position, so stamp only up to the deletion's real end
      uint32_t del_end = v.position - 1 + (v.nuc.size() - 1);
      if(del_end > i){
        is_del = true;
        deletion_span = del_end - i;
      }
      break;
    }

    for(auto &v : variant_records[i]){
      v.position_half_normal_upper = position_half_normal_upper;
    }
    if(is_del){
      for(uint32_t z=0; z < deletion_span && (i+z) < variant_records.size(); z++){
        if(variant_records[i+z].empty()){
          variant blank{};
          blank.position = i+z+1;
          blank.nuc = "-";
          blank.assigned_deletion = true;
          variant_records[i+z].push_back(blank);
          continue;
        }
        for(auto &v : variant_records[i+z]){
          v.assigned_deletion = true;
        }
      }
    }
  }
}

void consensus_sequence::write_consensus_to_file(std::string consensus_filename){
  std::ofstream file(consensus_filename, std::ios::app);
  std::string tmp = std::accumulate(sequence.begin(), sequence.end(), std::string(""));
  tmp.erase(std::remove(tmp.begin(), tmp.end(), '-'), tmp.end());
  file << ">"+seq_name << "\n";
  file << tmp << "\n";
  file.close();
}

void cluster_consensus(std::vector<variant> variants, \
                      std::string clustering_file, \
                      std::vector<double> solution){
  std::cerr << "calling consensus" << std::endl;
  if(variants.size() == 0){
    return;
  }

  std::string consensus_filename = clustering_file + ".fa";
 
  //find the largest position in the variants file
  uint32_t max_position = 0;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  }
  //initialize a consensus_sequence for each possible population
  std::vector<consensus_sequence> all_consensus_seqs;
  all_consensus_seqs.reserve(solution.size());
  for(uint32_t i=0; i < solution.size(); i++){
    all_consensus_seqs.emplace_back(max_position);
  }
  assign_variants_position(variants, all_consensus_seqs);
  //all_consensus_seqs[k] is indexed to match solution[k]/consensus_numbers, so this
  //processing loop stays in that original order; only the final write order changes.
  for(uint32_t i=0; i < all_consensus_seqs.size(); i++){
    all_consensus_seqs[i].set_seq_name(clustering_file + "_cluster_" + std::to_string(solution[i]));
    all_consensus_seqs[i].process_variant_assignments();
    all_consensus_seqs[i].get_consensus();
  }

  //write largest-abundance genome first: sort indices by solution value, descending
  std::vector<uint32_t> write_order(all_consensus_seqs.size());
  for(uint32_t i=0; i < write_order.size(); i++){
    write_order[i] = i;
  }
  std::sort(write_order.begin(), write_order.end(), [&](uint32_t a, uint32_t b){ return solution[a] > solution[b]; });

  std::ofstream(consensus_filename, std::ios::trunc); //start each run from an empty file, since write_consensus_to_file appends
  for(uint32_t idx : write_order){
    all_consensus_seqs[idx].write_consensus_to_file(consensus_filename);
  }
}
