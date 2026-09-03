#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <limits>
#include <set>
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/call_consensus_clustering.h"
#include "../src/solve_clustering.h"

void read_consensus(std::vector<std::pair<std::string, std::string>> &sequences, std::string filename){
  std::ifstream file(filename);
  std::string sequence;
  std::string name;
  std::string tmp;
  while (std::getline(file, tmp)) {
    if (tmp.find(">") != std::string::npos){
      if(sequence.size() > 0){
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),[](unsigned char c) { return std::toupper(c); });
        sequences.emplace_back(name, sequence);
      }
      name = tmp;
      sequence.clear();
      continue;
    }
    sequence += tmp;
  }
  if(sequence.size() > 0){
    sequences.emplace_back(name, sequence);
  }
}

int main() {
  std::string prefix = "/tmp/consensus";
  int num_tests = 4;
  int success = 0;

  uint32_t min_depth = 5;
  uint8_t min_qual = 20;
  double default_threshold = 0.5;
  double invariant_threshold = 0.97;

  uint32_t n = 4;

  //TEST 1 - manually currated data
  std::string var_filename = "../data/version_bump_tests/vbump_consensus_var.txt";
  std::string consensus_filename = "../data/version_bump_tests/vbump_consensus_ivar.fa";

  std::vector<double> solution;
  std::vector<double> means;
  std::vector<variant> variants = gmm_model(var_filename, prefix, min_depth, min_qual, solution, means, default_threshold, n, invariant_threshold, 1e-3, 1e-2, 0.05);

  // TEST 2 - verify the dominant deletion -GGA at POS=486 (stored as position=487)
  // is not excluded as an overlapped deletion and is assigned to at least one
  // consensus with no duplicate consensus_numbers.
  {
    bool found = false;
    for (const auto& v : variants) {
      if (v.position == 487 && v.nuc == "-GGA") {
        found = true;
        std::set<uint32_t> cn_set(v.consensus_numbers.begin(), v.consensus_numbers.end());
        bool no_dups = cn_set.size() == v.consensus_numbers.size();
        std::cerr << "[GGA check] overlapped_deletion=" << v.overlapped_deletion
                  << " consensus_numbers=[ ";
        for (auto cn : v.consensus_numbers) std::cerr << cn << " ";
        std::cerr << "] cluster=" << v.cluster_assigned
                  << " gapped_freq=" << v.gapped_freq << std::endl;
        if (!v.overlapped_deletion && !v.consensus_numbers.empty() && no_dups) {
          success++;
        } else {
          std::cerr << "-GGA consensus_numbers check failed" << std::endl;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "-GGA deletion at position 487 not found in variants" << std::endl;
    }
  }

  // TEST 3 - verify the dominant insertion +TT at POS=1553 (gapped_freq=0.9) is
  // assigned to the major/dominant consensus genome (the one with the largest
  // accepted solution fraction).
  {
    bool found = false;
    for (const auto& v : variants) {
      if (v.position == 1553 && v.nuc == "+TT") {
        found = true;
        uint32_t major_idx = std::distance(solution.begin(), std::max_element(solution.begin(), solution.end()));
        bool in_major = std::find(v.consensus_numbers.begin(), v.consensus_numbers.end(), major_idx) != v.consensus_numbers.end();
        std::cerr << "[+TT check] gapped_freq=" << v.gapped_freq << " cluster=" << v.cluster_assigned
                  << " major_idx=" << major_idx << " consensus_numbers=[ ";
        for (auto cn : v.consensus_numbers) std::cerr << cn << " ";
        std::cerr << "]" << std::endl;
        if (in_major) {
          success++;
        } else {
          std::cerr << "+TT insertion at position 1553 not assigned to major genome" << std::endl;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "+TT insertion at position 1553 not found in variants" << std::endl;
    }
  }

  // TEST 4 - a genome that has a confident, single-nucleotide "confirmed" call
  // at a position must still be forced to N if a different explanation for
  // that same position also implicates it as merely ambiguous (the classic
  // 10%/30% shared-mutation collision landing on the same peak as a private
  // 40% mutation). Tested directly against assign_variants_position/
  // get_consensus so it's independent of how solve_clustering.cpp derived
  // consensus_numbers/ambiguous_numbers - that combinatorics is covered
  // separately in test_solve_clustering.cpp.
  {
    variant v_confirmed;
    v_confirmed.position = 5;
    v_confirmed.nuc = "A";
    v_confirmed.consensus_numbers = {3}; //genome 3 looks confidently assigned here

    variant v_ambiguous;
    v_ambiguous.position = 5;
    v_ambiguous.nuc = "T";
    v_ambiguous.ambiguous_numbers = {3}; //but the same genome is also implicated by a colliding explanation

    variant v_private;
    v_private.position = 8;
    v_private.nuc = "G";
    v_private.consensus_numbers = {1}; //unrelated genome, unambiguous - should be unaffected

    std::vector<variant> ambiguity_variants = {v_confirmed, v_ambiguous, v_private};

    std::vector<consensus_sequence> ambiguity_consensus_seqs;
    for(uint32_t i=0; i < 4; i++){
      ambiguity_consensus_seqs.emplace_back(10);
    }
    assign_variants_position(ambiguity_variants, ambiguity_consensus_seqs);
    for(uint32_t i=0; i < ambiguity_consensus_seqs.size(); i++){
      ambiguity_consensus_seqs[i].process_variant_assignments();
      ambiguity_consensus_seqs[i].get_consensus(i);
    }

    bool forced_ambiguous = ambiguity_consensus_seqs[3].get_base(5) == "N";
    bool private_still_confident = ambiguity_consensus_seqs[1].get_base(8) == "G";

    if(forced_ambiguous && private_still_confident){
      success++;
    } else {
      std::cerr << "TEST 4 failed, genome 3 pos 5 = " << ambiguity_consensus_seqs[3].get_base(5)
                << " genome 1 pos 8 = " << ambiguity_consensus_seqs[1].get_base(8) << std::endl;
    }
  }

  if(!variants.empty()) {
    cluster_consensus(variants, prefix, default_threshold, min_depth, min_qual, solution, means);
  }

  uint32_t min_position = std::numeric_limits<uint32_t>::max();
  for (const auto& v : variants) {
    if (v.total_depth > 0 && v.position < min_position)
      min_position = v.position;
  }

  std::vector<std::pair<std::string, std::string>> gt_sequences;
  read_consensus(gt_sequences, consensus_filename);
  std::string exp_sequence;
  std::vector<std::pair<std::string, std::string>> exp_sequences;
  read_consensus(exp_sequences, prefix+".fa");

  // Build a map from trimmed-string index to genomic position accounting for deletions.
  // Deletions remove characters from the consensus string, so index i in the trimmed
  // string does NOT equal min_position + i when deletions fall within the range.
  std::set<uint32_t> deleted_positions;
  for (const auto& v : variants) {
    if (v.nuc.find('-') == std::string::npos) continue;
    std::string nuc = v.nuc;
    nuc.erase(std::remove(nuc.begin(), nuc.end(), '-'), nuc.end());
    for (uint32_t z = 0; z < nuc.size(); z++) {
      deleted_positions.insert(v.position + z);
    }
  }
  // Build index->genomic mapping: production output (exp_sequences) is no longer
  // trimmed to min_position, so it starts at genomic position 1. Walk from there,
  // skip deleted ones, and assign each surviving position to the next string index.
  std::vector<uint32_t> index_to_genomic;
  {
    uint32_t max_pos = 0;
    for (const auto& v : variants) if (v.position > max_pos) max_pos = v.position;
    for (uint32_t p = 1; p <= max_pos + 10; p++) {
      if (deleted_positions.count(p) == 0) {
        index_to_genomic.push_back(p);
      }
    }
  }
  // gt_sequences is trimmed and starts at min_position; find where that lands in
  // the untrimmed exp_sequences index mapping so the two can be compared aligned.
  uint32_t exp_start = 0;
  while (exp_start < index_to_genomic.size() && index_to_genomic[exp_start] < min_position) {
    exp_start++;
  }

  bool correct = true;
  for (auto itgt = gt_sequences.begin(), itexp = exp_sequences.begin(); itgt != gt_sequences.end() && itexp != exp_sequences.end(); ++itgt, ++itexp) {
    uint32_t exp_len = itexp->second.size() > exp_start ? (uint32_t)itexp->second.size() - exp_start : 0;
    if(exp_len != itgt->second.size()) {
      correct = false;
      std::cerr << "not same size gt=" << itgt->second.size() << " exp=" << exp_len << std::endl;
    }
    for(uint32_t i=0; i < exp_len && i < itgt->second.size(); i++){
      char a = itgt->second[i];
      char b = itexp->second[exp_start + i];

      if(a != b){
        uint32_t genomic_pos = (exp_start+i < index_to_genomic.size()) ? index_to_genomic[exp_start+i] : (min_position + i);
        std::cerr << "gt " << a << " exp " << b << " consensus_pos=" << i+1 << " genomic_pos=" << genomic_pos << std::endl;
        exit(0);
        int ctx = 5;
        std::cerr << "  gt  context: ";
        for (int c = -ctx; c <= ctx; c++) {
          int idx = (int)i + c;
          if (idx < 0 || idx >= (int)itgt->second.size()) std::cerr << '.';
          else std::cerr << (c == 0 ? '[' : ' ') << itgt->second[idx] << (c == 0 ? ']' : ' ');
        }
        std::cerr << std::endl;
        std::cerr << "  exp context: ";
        for (int c = -ctx; c <= ctx; c++) {
          int idx = (int)(exp_start + i) + c;
          if (idx < 0 || idx >= (int)itexp->second.size()) std::cerr << '.';
          else std::cerr << (c == 0 ? '[' : ' ') << itexp->second[idx] << (c == 0 ? ']' : ' ');
        }
        std::cerr << std::endl;
        for (const auto& v : variants) {
          if (v.position == genomic_pos) {
            std::cerr << "  nuc=" << v.nuc
                      << " gapped_freq=" << v.gapped_freq
                      << " gapped_depth=" << v.gapped_depth
                      << " qual=" << v.qual
                      << " cluster=" << v.cluster_assigned
                      << " half_normal_upper=" << v.half_normal_upper
                      << " half_normal_lower=" << v.half_normal_lower
                      << " position_masked=" << v.position_masked
                      << " amplicon_masked=" << v.amplicon_masked
                      << " depth_flag=" << v.depth_flag
                      << " qual_flag=" << v.qual_flag
                      << " outside_freq_range=" << v.outside_freq_range
                      << " consensus_numbers=[ ";
            for (auto cn : v.consensus_numbers) std::cerr << cn << " ";
            std::cerr << "]" << std::endl;
          }
        }
        correct = false;
      }
    }
  }
  if(correct) success++;

  std::cerr << "num tests " << num_tests << " success " << success << std::endl;
  return (num_tests == success) ? 0 : -1;
}
