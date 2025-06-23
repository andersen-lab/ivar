#include "get_masked_amplicons.h"

int get_primers_with_mismatches(std::string bed, std::string vpath,
                                std::string out, std::string primer_pair_file) {
  std::vector<primer> primers = populate_from_file(bed);
  std::vector<primer> mismatches_primers;
  std::vector<primer> tmp;
  if (primers.size() == 0) {
    return 0;
  }
  populate_pair_indices(primers, primer_pair_file);
  std::string line, cell;
  std::ifstream fin(vpath.c_str());
  out += ".txt";
  std::ofstream fout(out.c_str());
  unsigned int ctr, pos;
  std::stringstream line_stream;
  while (std::getline(fin, line)) {
    line_stream << line;
    ctr = 0;
    pos = 0;
    while (std::getline(line_stream, cell, '\t')) {
      switch (ctr) {
        case 1:
          if (cell != "POS") pos = stoi(cell);
          break;
        default:
          break;
      }
      ctr++;
    }
    if (pos == 0) {
      line_stream.clear();
      continue;
    }
    pos--;  // 1 based to 0 based
    tmp = get_primers(primers, pos);
    // mismatches_primers.insert(mismatches_primers.end(), tmp.begin(),
    // tmp.end());
    std::vector<primer>::iterator tmp_it;
    for (std::vector<primer>::iterator it = tmp.begin(); it != tmp.end();
         ++it) {
      tmp_it =
          std::find(mismatches_primers.begin(), mismatches_primers.end(), *it);
      if (tmp_it == mismatches_primers.end()) {
        std::cout << it->get_name() << std::endl;
        mismatches_primers.push_back(*it);
      }
      // Look for primer pair
      if (it->get_pair_indice() != -1) {
        tmp_it = std::find(mismatches_primers.begin(), mismatches_primers.end(),
                           primers.at(it->get_pair_indice()));
        if (tmp_it == mismatches_primers.end()) {
          mismatches_primers.push_back(primers.at(it->get_pair_indice()));
        }
      }
    }
    line_stream.clear();
    tmp.clear();
  }
  for (std::vector<primer>::iterator it = mismatches_primers.begin();
       it != mismatches_primers.end(); ++it) {
    fout << it->get_name();
    std::cout << it->get_name();
    if (it != mismatches_primers.end() - 1) {
      fout << "\t";
      std::cout << "\t";
    }
  }
  std::cout << std::endl;
  return 0;
}
