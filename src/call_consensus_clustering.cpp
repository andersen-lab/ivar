#include "estimate_error.h"
#include "call_consensus_clustering.h"
#include "gmm.h"
#include "saga.h"
#include "ref_seq.h"
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>

std::string trim_leading_ambiguities(std::string sequence, uint32_t min_position){
  std::string result = sequence.substr(min_position-1);
  return(result);
}

void call_majority_consensus(std::vector<variant> variants, uint32_t max_position, std::string clustering_file, double default_threshold){
  //if we can't find a solution simply take the majority variant per position
  std::vector<std::string> nucs;
  std::vector<double> freqs;
  std::vector<std::string> tmp(max_position, "N");
  for(uint32_t i=1; i <= max_position; i++){
    freqs.clear();
    nucs.clear();
    for(uint32_t j=0; j < variants.size(); j++){
      if(variants[j].position == i){
        nucs.push_back(variants[j].nuc);
        freqs.push_back(variants[j].freq);
      }
    }
    if(freqs.size() == 0) continue;
    uint32_t index = std::distance(freqs.begin(), std::max_element(freqs.begin(), freqs.end()));
    if(freqs[index] >= (double)default_threshold){
      tmp[i-1] = nucs[index];
    }
  }
  std::string consensus_string = std::accumulate(tmp.begin(), tmp.end(), std::string(""));
  //write the consensus to file
  std::string consensus_filename = clustering_file + ".fa";
  std::ofstream file(consensus_filename);
  std::string name = ">"+clustering_file+"_"+std::to_string(default_threshold)+"_threshold";
  file << name << "\n";
  file << consensus_string << "\n";
  file.close();
}

void cluster_consensus(std::vector<variant> variants, std::string clustering_file, double default_threshold, uint32_t min_depth, uint8_t min_qual, std::vector<double> solution, std::vector<double> means, std::string ref){
  //TODO call majority
  if(variants.size() == 0){
    std::cerr << "haven't solved this yet" << std::endl;
  }

  std::cerr << "calling consensus" << std::endl;
  //parse reference sequence
  ref_antd refantd(ref, "");

  double max_mean=0;
  set_freq_range_flags(variants, 0, 1);
  double error_rate = cluster_error(variants, min_qual, min_depth);
  double freq_lower_bound = 1-error_rate+0.001;
  double freq_upper_bound = error_rate-0.001;
  set_freq_range_flags(variants, freq_lower_bound, freq_upper_bound);
  //find the largest position in the variants file
  uint32_t max_position = 0;
  uint32_t min_position = 4294967295U;;
  for(auto x : variants){
    if(x.position > max_position){
      max_position = x.position;
    }
    if(x.position < min_position && x.total_depth > 0){
        min_position = x.position;
    }
  }
  bool print = false;
  //initialize sequences for all possible populations
  std::vector<std::vector<std::string>> all_consensus_seqs;
  for(uint32_t i=0; i < means.size(); i++){
    std::vector<std::string> tmp(max_position, "N");
    all_consensus_seqs.push_back(tmp);
  }

  for(auto m  : means){
    std::cerr << m << " ";
  }
  //order varaints by position
  std::sort(variants.begin(), variants.end(), [](const variant& a, const variant& b) {return a.position < b.position;});
  std::vector<uint32_t> last_adjustment(all_consensus_seqs.size(), 0);

  //track deletions over time
  std::vector<std::vector<uint32_t>> deletions(means.size());

  //iterate all variants and determine
  for(uint32_t i = 0; i < variants.size(); i++){
    //TESTLINES
    if(variants[i].position == 489){
      print = true;
      std::cerr << "\ntop freq " << variants[i].freq << " " << variants[i].nuc << " cluster " << variants[i].cluster_assigned << " gapped freq " << variants[i].gapped_freq << std::endl;
      std::cerr << "vague assignment " << variants[i].vague_assignment << " depth flag " << variants[i].depth_flag << std::endl;
      std::cerr << "amplicon masked " << variants[i].amplicon_masked << " amp flux pos " << variants[i].amplicon_flux << std::endl;
    }else{
      print = false;
    }
    double freq = variants[i].gapped_freq;
    double qual = variants[i].qual;
    uint32_t depth = variants[i].gapped_depth;
    //depth, quality, and low frequency bypass
    if(freq < freq_lower_bound || qual < (double)min_qual || depth < min_depth){
      if(print) std::cerr << "min qual, freq, or depth issue " << qual << " " << freq << " " << depth << " flb " << freq_lower_bound << " mq " << (double)min_qual << " md " << min_depth << std::endl;
      continue;
    }

    //if this amplicon is experiencing fluctuation across amplicons, call ambiguity
    if(variants[i].amplicon_masked && variants[i].freq < freq_upper_bound){
      if(print){
        std::cerr << "amplicon is experiencing fluctuation" << std::endl;
      }
      continue;
    }

    //this variant position experiences fluctuation across amplicons
    if(variants[i].amplicon_flux){
      if(print){
        std::cerr << "amplicon in flux" << std::endl;
      }
      continue;
    }
     uint32_t position = variants[i].position;
     if(variants[i].vague_assignment && variants[i].freq < freq_upper_bound && variants[i].freq < max_mean){
       if(print){
          std::cerr << "d" << std::endl;
          for(auto a : variants[i].probabilities){
            std::cerr << a << " ";
          }
          std::cerr << "\n";
       }
       continue;
     }
     bool del = variants[i].nuc.find('-') != std::string::npos;
     //handle all the cases where you never assigned anything, assign to all if it's over the upper bound
     if(variants[i].cluster_assigned == -1){
      if(variants[i].gapped_freq < freq_upper_bound) continue;
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
  std::cerr << "printing deletions" << std::endl;
  for(auto con : deletions){
    for(auto del : con){
      std::cerr << del << " ";
    }
    std::cerr << "\n";
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
    std::cerr << "mean " << tmp_mean << std::endl;
    auto it = std::find(solution.begin(), solution.end(), tmp_mean);
    if(it == solution.end()){
      continue;
    }
    for(uint32_t j=485; j < 490; j++){
      std::cerr << j << " " << sorted_strings[i][j] << " ";
    }
    std::cerr << "\n";
    std::string trimmed_sequence = trim_leading_ambiguities(sorted_strings[i], min_position);
    file << ">"+clustering_file+"_cluster_"+ std::to_string(tmp_mean) << "\n";
    file << trimmed_sequence << "\n";
  }
  file.close();
}
