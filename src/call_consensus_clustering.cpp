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

void call_majority_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold){
  std::cerr << "in majority consensus call" << std::endl;
  uint32_t max_position = 0;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
  }
  consensus_sequence cs(max_position); 
  std::vector<consensus_sequence> all_consensus_seqs = {cs};
  assign_variants_position(variants, all_consensus_seqs);
  all_consensus_seqs[0].set_seq_name(clustering_file + "_" + std::to_string(default_threshold) + "_threshold");
  all_consensus_seqs[0].process_variant_assignments();
  all_consensus_seqs[0].get_majority_consensus(default_threshold);
  std::string consensus_filename = clustering_file + "_threshold.fa";
  std::ofstream(consensus_filename, std::ios::trunc); //start each run from an empty file, since write_consensus_to_file appends
  all_consensus_seqs[0].write_consensus_to_file(consensus_filename);
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

void consensus_sequence::get_consensus(uint32_t n){
  uint32_t test_position = 4184-1;
  uint32_t deletion_span; //track deletion spans
  for(uint32_t i=0; i < variant_records.size(); i++){
    if(position_ambiguous[i]){
      continue; //multiple explanations for this peak - leave as N
    }
    if(variant_records[i].size() == 0){
      continue;
    }
    std::vector<string> nucs; //all the possible nucs here?
    bool deletion_added = false;
    std::string insertion;
    for(uint32_t j=0; j < variant_records[i].size(); j++){
      /*if(i == test_position){
        std::cerr << "\nconsensus " << n << " position " << i+1 << " nuc " << variant_records[i][j].nuc << " freq " << variant_records[i][j].gapped_freq << std::endl;
        std::cerr << "upper half normal " << variant_records[i][j].half_normal_upper << " lower half normal " << variant_records[i][j].half_normal_lower << std::endl;
        std::cerr << "assigned deletion " << variant_records[i][j].assigned_deletion << std::endl;
        for(auto cc : variant_records[i][j].consensus_numbers){
          std::cerr << cc << " ";
        }
        std::cerr << std::endl;
      }*/

      bool has_insertion = variant_records[i][j].nuc.find('+') != std::string::npos;
        if(has_insertion){
          insertion = variant_records[i][j].nuc;
          insertion.erase(std::remove(insertion.begin(), insertion.end(), '+'), insertion.end());
      }

      //if we have one assigned to the upper half normal, we go with that and ignore the rest

      if(variant_records[i][j].half_normal_upper){
        nucs.push_back(variant_records[i][j].nuc);
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
      sequence[i] = std::string(1, combined);
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

    for(uint32_t j=0; j < variant_records[i].size(); j++){
      if(variant_records[i][j].qual_flag || variant_records[i][j].depth_flag){
        continue;
      }

      bool has_insertion = variant_records[i][j].nuc.find('+') != std::string::npos;
      if(has_insertion){
        insertion = variant_records[i][j].nuc;
        insertion.erase(std::remove(insertion.begin(), insertion.end(), '+'), insertion.end());
        continue;
      }

      if(variant_records[i][j].assigned_deletion){
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
    bool position_half_normal_upper = false;
    bool is_del = false;
    uint32_t deletion_span = 0;

    for(uint32_t j=0; j < variant_records[i].size(); j++){
      if(variant_records[i][j].half_normal_upper){
        position_half_normal_upper = true;
        break;
      }
      if(variant_records[i][j].nuc.find('-') != std::string::npos){
        is_del = true;
        deletion_span = variant_records[i][j].nuc.size()-1;
        break;
      }
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
                      double default_threshold, \
                      uint32_t min_depth, \
                      uint8_t min_qual, \
                      std::vector<double> solution, \
                      std::vector<double> means){
  std::cerr << "calling consensus" << std::endl;
  if(variants.size() == 0){
    return;
  }

  std::string consensus_filename = clustering_file + ".fa";
 
  //find the largest position in the variants file
  uint32_t max_position = 0;
  uint32_t min_position = 4294967295U;;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
    if(x.position < min_position && x.gapped_depth > 0){
        min_position = x.position;
    }
  }
  bool print = false;
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
    all_consensus_seqs[i].get_consensus(i);
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
  /*std::vector<uint32_t> last_adjustment(all_consensus_seqs.size(), 0);

  //track deletions over time
  std::vector<std::vector<uint32_t>> deletions(means.size());
  //iterate all variants and determine
  for(uint32_t i = 0; i < variants.size(); i++){
    //TESTLINES
    if(variants[i].position == 21770){
      print = true;
      std::cerr << "\ntop freq " << variants[i].freq << " " << variants[i].nuc << " cluster " << variants[i].cluster_assigned << " gapped freq " << variants[i].gapped_freq << " depth " << variants[i].total_depth << std::endl;
      std::cerr << "vague assignment " << variants[i].vague_assignment << " depth flag " << variants[i].depth_flag << std::endl;
      std::cerr << "amplicon masked " << variants[i].amplicon_masked << " position masked " << variants[i].position_masked << std::endl;
      std::cerr << "cluster assigned " << variants[i].cluster_assigned << std::endl;
      for(auto cc : variants[i].consensus_numbers){
        std::cerr << "consensus number " << cc << std::endl;
      }
      std::cerr << "half normal upper " << variants[i].half_normal_upper << " half normal lower " << variants[i].half_normal_lower << std::endl;
    }else{
      print = false;
    }

    if(variants[i].position_conflict){
      if(print) std::cerr << "position conflict " << variants[i].position << std::endl;
      continue;
    }
    if(variants[i].overlapped_deletion){
      continue;
    }
    if(variants[i].position_masked && !variants[i].half_normal_upper){
      if(print) std::cerr << "position masked " << variants[i].position_masked << std::endl;
      continue;
    }
    double freq = variants[i].gapped_freq;
    double qual = variants[i].qual;
    uint32_t depth = variants[i].gapped_depth;
    //depth, quality, and low frequency bypass
    if(variants[i].half_normal_lower || qual < (double)min_qual || depth < min_depth || variants[i].vague_assignment){
      if(print) {
        std::cerr << "half normal lower " << variants[i].half_normal_lower << std::endl;
        std::cerr << "min qual or depth issue " << qual << " "  << depth <<  " mq " << (double)min_qual << " md " << min_depth << std::endl;
        std::cerr << "vague assignment " << variants[i].vague_assignment << std::endl;
      }
      continue;
    }

    //if this amplicon is experiencing fluctuation across amplicons, call ambiguity
    if(variants[i].amplicon_masked && !variants[i].half_normal_upper){
      if(print){
        std::cerr << "amplicon is experiencing fluctuation" << std::endl;
      }
      continue;
    }
    }

    uint32_t position = variants[i].position;
    if(variants[i].vague_assignment && !variants[i].half_normal_upper){
       if(print){
          for(auto a : variants[i].probabilities){
            std::cerr << a << " ";
          }
          std::cerr << "\n";
       }
       continue;
    }

    bool del = variants[i].nuc.find('-') != std::string::npos;

    //handle all the cases where you never assigned anything, assign to all if it's over the upper bound
     if(variants[i].half_normal_upper){
      if(print) std::cerr << "not assigned anything" << std::endl;
      for(uint32_t j=0; j < all_consensus_seqs.size(); j++){
        uint32_t adjusted_pos = position-1;
        if(variants[i].nuc.find('+') != std::string::npos){
          std::string nuc = variants[i].nuc;
          nuc.erase(std::remove(nuc.begin(), nuc.end(), '+'), nuc.end());
          if(last_adjustment[j] == position){
            all_consensus_seqs[j][position-1] += nuc;
          } else {
            all_consensus_seqs[j][position-1] = nuc;
          }
          last_adjustment[j] = position;
        } else if (variants[i].position == last_adjustment[j] && !del){
          all_consensus_seqs[j][adjusted_pos].insert(0, variants[i].nuc);
        } else {
          if(!del){
            if(print) std::cerr << "here " << j << " " << adjusted_pos << " " << variants[i].nuc << std::endl;
            all_consensus_seqs[j][adjusted_pos] = variants[i].nuc;
            last_adjustment[j] = position;
          } else {
            std::string nuc = variants[i].nuc;
            nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
            for(uint32_t z=0; z < nuc.size(); z++){
              all_consensus_seqs[j][position-1+z] = "-";
              deletions[j].push_back(position-1+z);
            }
          }
        }
      }
      continue;
    }

    for(uint32_t j=0; j < variants[i].consensus_numbers.size(); j++){
        uint32_t k = variants[i].consensus_numbers[j];
        bool found_del = std::find(deletions[k].begin(), deletions[k].end(), variants[i].position) != deletions[k].end();
        if(found_del) continue; //already assigned a deletion to this position
        if(variants[i].nuc.find('+') != std::string::npos){
          std::string nuc = variants[i].nuc;
          nuc.erase(std::remove(nuc.begin(), nuc.end(), '+'), nuc.end());
          if(last_adjustment[k] == position){
            all_consensus_seqs[k][position-1] += nuc;
          } else {
            all_consensus_seqs[k][position-1] = nuc;
          }
          last_adjustment[k] = position;
        } else if (variants[i].position == last_adjustment[k] && !del){
          all_consensus_seqs[k][position-1].insert(0, variants[i].nuc);
        } else {
          if(!del){
            all_consensus_seqs[k][position-1] = variants[i].nuc;
            last_adjustment[k] = position;
          } else {
            std::string nuc = variants[i].nuc;
            nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
            for(uint32_t z=0; z < nuc.size(); z++){
              all_consensus_seqs[k][position-1+z] = "-";
              deletions[k].push_back(position-1+z);
            }
          }
      }
    }
  }

  std::vector<std::string> all_sequences;
  for(uint32_t i=0; i < all_consensus_seqs.size(); i++){
    std::string tmp = std::accumulate(all_consensus_seqs[i].begin(), all_consensus_seqs[i].end(), std::string(""));
    tmp.erase(std::remove(tmp.begin(), tmp.end(), '-'), tmp.end());
    all_sequences.push_back(tmp);
  }
  //write the consensus string to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);

  std::vector<uint32_t> indices(all_sequences.size());
  for (uint32_t i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }

  // Sort indices based on double values (descending)
  std::sort(indices.begin(), indices.end(), [&](uint32_t i, uint32_t j) {return means[i] > means[j];});

  // Apply sorted order
  std::vector<std::string> sorted_strings;
  std::vector<double> sorted_values;
  for (auto i : indices) {
    sorted_strings.push_back(all_sequences[i]);
    sorted_values.push_back(means[i]);
  }

  for(uint32_t i=0; i < sorted_strings.size(); i++){
    double tmp_mean = sorted_values[i];
    auto it = std::find(solution.begin(), solution.end(), tmp_mean);
    if(it == solution.end()){
      continue;
    }

    std::string next_trimmed_sequence = trim_leading_ambiguities(sorted_strings[i], min_position);
    file << ">"+clustering_file+"_cluster_"+ std::to_string(tmp_mean) << "\n";
    file << next_trimmed_sequence << "\n";
  }
  file.close();*/
}
