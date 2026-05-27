#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../src/site_state.h"
#include "../src/variant_caller.h"

static site_state make_nuc(char base, uint32_t pos, uint8_t qual = 30) {
  return site_state(base, qual, pos, NUCLEOTIDE);
}

static bool file_exists(const std::string &path) {
  return std::ifstream(path).good();
}

static std::vector<std::string> split_tab(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) fields.push_back(field);
  return fields;
}

int main() {
  int success = 0;
  const std::string out_prefix = "/tmp/test_write_codon_aa";
  const std::string out_prefix_nogff = "/tmp/test_write_codon_aa_nogff";

  // End-to-end: 3 TGA (stop), 7 AGC (ref).
  variant_caller vc(20, "../../data/db/test_ref.fa", "../../data/test_discontinuous_cds.gff");
  vc.initialize_region("test");

  for (int r = 0; r < 3; ++r) {
    std::vector<site_state> states;
    states.push_back(make_nuc('T', 10));
    states.push_back(make_nuc('G', 11));
    states.push_back(make_nuc('A', 12));
    vc.add_variants(states);
  }
  for (int r = 0; r < 7; ++r) {
    std::vector<site_state> states;
    states.push_back(make_nuc('A', 10));
    states.push_back(make_nuc('G', 11));
    states.push_back(make_nuc('C', 12));
    vc.add_variants(states);
  }

  vc.write_codon_to_file(out_prefix, "test");
  vc.write_aa_to_file(out_prefix, "test");

  if (file_exists(out_prefix + ".codons.txt"))
    success += 1;
  else
    std::cout << "Codon file not created" << std::endl;

  if (file_exists(out_prefix + ".aa.txt"))
    success += 1;
  else
    std::cout << "AA file not created" << std::endl;

  // Codon file header + TGA row content.
  std::ifstream codon_file(out_prefix + ".codons.txt");
  std::string header;
  std::getline(codon_file, header);
  const std::string &expected_codon_header = variant_caller::CODON_FILE_HEADER;
  if (header == expected_codon_header)
    success += 1;
  else
    std::cout << "Codon header mismatch: '" << header << "'" << std::endl;

  bool found_tga = false;
  std::string line;
  while (std::getline(codon_file, line)) {
    std::vector<std::string> f = split_tab(line);
    if (f.size() < 10)
      continue;
    if (f[2] == "AGC" && f[3] == "TGA") {
      found_tga = true;
      if (f[1] == "1")
        success += 1;
      else
        std::cout << "Codon POS_CODON mismatch: '" << f[1] << "'" << std::endl;
      if (f[7] == "10,12")
        success += 1;
      else
        std::cout << "Codon POS list mismatch: '" << f[7] << "'" << std::endl;
      if (f[8] == "A,C")
        success += 1;
      else
        std::cout << "Codon REF list mismatch: '" << f[8] << "'" << std::endl;
      if (f[9] == "T,A")
        success += 1;
      else
        std::cout << "Codon ALT list mismatch: '" << f[9] << "'" << std::endl;
      if (f[6] == "0.3")
        success += 1;
      else
        std::cout << "Codon ALT_FREQ_CODON mismatch: '" << f[6] << "'" << std::endl;
      break;
    }
  }
  if (found_tga)
    success += 1;
  else
    std::cout << "TGA row not found in codon file" << std::endl;

  // AA file header + stop AA row content.
  std::ifstream aa_file(out_prefix + ".aa.txt");
  std::getline(aa_file, header);
  const std::string &expected_aa_header = variant_caller::AA_FILE_HEADER;
  if (header == expected_aa_header)
    success += 1;
  else
    std::cout << "AA header mismatch: '" << header << "'" << std::endl;

  bool found_stop = false;
  while (std::getline(aa_file, line)) {
    std::vector<std::string> f = split_tab(line);
    if (f.size() < 9)
      continue;
    if (f[3] == "*") {
      found_stop = true;
      if (f[1] == "1")
        success += 1;
      else
        std::cout << "AA POS_AA mismatch: '" << f[1] << "'" << std::endl;
      if (f[8] == "TGA")
        success += 1;
      else
        std::cout << "AA ALT_CODON mismatch: '" << f[8] << "'" << std::endl;
      if (f[6] == "0.3")
        success += 1;
      else
        std::cout << "AA ALT_FREQ_AA mismatch: '" << f[6] << "'" << std::endl;
      break;
    }
  }
  if (found_stop) ++success;
  else std::cerr << "Stop AA row not found" << std::endl;

  // No GFF: write should produce no files
  variant_caller vc2(20, "../../data/db/test_ref.fa");
  vc2.initialize_region("test");

  std::vector<site_state> states;
  states.push_back(make_nuc('T', 10));
  states.push_back(make_nuc('G', 11));
  states.push_back(make_nuc('A', 12));
  vc2.add_variants(states);

  vc2.write_codon_to_file(out_prefix_nogff, "test");
  vc2.write_aa_to_file(out_prefix_nogff, "test");

  if (!file_exists(out_prefix_nogff + ".codons.txt"))
    success += 1;
  else
    std::cerr << "Codon file should NOT exist without GFF" << std::endl;
  if (!file_exists(out_prefix_nogff + ".aa.txt"))
    success += 1;
  else
    std::cerr << "AA file should NOT exist without GFF" << std::endl;

  return success == 16 ? 0 : 1;
}
