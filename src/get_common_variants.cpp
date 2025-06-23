#include "get_common_variants.h"

const int NUM_FIELDS = 19;
const std::string fields[NUM_FIELDS] = {
    "REGION",    "POS",      "REF",       "ALT",    "REF_DP",
    "REF_RV",    "REF_QUAL", "ALT_DP",    "ALT_RV", "ALT_QUAL",
    "ALT_FREQ",  "TOTAL_DP", "PVAL",      "PASS",   "GFF_FEATURE",
    "REF_CODON", "REF_AA",   "ALT_CODON", "ALT_AA"};

const std::string na_tab_delimited_str =
    "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";

int read_variant_file(
    std::ifstream &fin, unsigned int file_number,
    std::map<std::string, unsigned int> &counts,
    std::map<std::string, std::string> &file_tab_delimited_str) {
  unsigned int ctr;
  std::string line, cell, tab_delimited_key, tab_delimited_val;
  std::stringstream line_stream;
  std::string region, ref, alt;
  // Make sure format of header matches
  std::getline(fin, line);
  line_stream << line;
  ctr = 0;
  while (std::getline(line_stream, cell, '\t')) {
    if (cell.compare(fields[ctr]) != 0) {
      return -1;
    }
    ctr++;
  }
  line_stream.clear();
  while (std::getline(fin, line)) {
    line_stream << line;
    ctr = 0;
    tab_delimited_key = "";
    tab_delimited_val = "";
    while (std::getline(line_stream, cell, '\t')) {
      switch (ctr) {
        case 0:   // REGION
        case 1:   // POS
        case 2:   // REF
        case 3:   // ALT
        case 14:  // GFF_FEATURE
        case 15:  // REF_CODON
        case 16:  // REF_AA
        case 17:  // ALT_CODON
        case 18:  // ALT_AA
          tab_delimited_key += cell + "\t";
          break;
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        default:
          tab_delimited_val += cell;
          if (ctr < 13) {
            tab_delimited_val += "\t";
          }
          break;
      }
      ctr++;
    }
    if (counts.find(tab_delimited_key) == counts.end()) {
      counts[tab_delimited_key] = 1;
    } else {
      counts[tab_delimited_key] += 1;
    }
    file_tab_delimited_str[tab_delimited_key + std::to_string(file_number)] =
        tab_delimited_val;
    line_stream.clear();
  }
  return 0;
}

int common_variants(std::string out, double min_threshold, char *files[],
                    unsigned int nfiles) {
  out += ".tsv";
  std::ofstream fout(out.c_str());
  unsigned int i, j;
  std::ifstream fin;
  std::map<std::string, unsigned int> counts =
      std::map<std::string, unsigned int>();
  std::map<std::string, std::string> file_tab_delimited_str =
      std::map<std::string, std::string>();
  for (i = 0; i < nfiles; ++i) {
    fin.open(files[i]);
    if (read_variant_file(fin, i, counts, file_tab_delimited_str) != 0) {
      std::cout << "Header format of " << files[i] << " did not match!";
      std::cout
          << " Please use files generated using \"ivar variants\" command."
          << std::endl;
    }
    fin.close();
  }
  std::map<std::string, unsigned int>::iterator it = counts.begin();
  // Write header
  for (i = 0; i < 4; ++i) {
    fout << fields[i] << "\t";
  }
  for (i = 14; i < 19; ++i) {
    fout << fields[i] << "\t";
  }
  for (i = 0; i < nfiles; ++i) {
    for (j = 4; j < 14; ++j) {
      fout << fields[j] << "_" << files[i] << "\t";
    }
  }
  fout << "\n";
  // Write rows
  while (it != counts.end()) {
    if (((float)it->second) / (float)nfiles >=
        min_threshold) {  // Check if variant occurs in more 'min_threshold'
                          // fraction of files
      fout << it->first;
      for (j = 0; j < nfiles; ++j) {
        if (file_tab_delimited_str.find(it->first + std::to_string(j)) ==
            file_tab_delimited_str.end()) {
          fout << na_tab_delimited_str;
        } else {
          fout << file_tab_delimited_str[it->first + std::to_string(j)];
        }
        if (j < nfiles - 1) fout << "\t";
      }
      fout << "\n";
    }
    it++;
  }
  return 0;
}
