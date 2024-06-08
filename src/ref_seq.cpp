#include "ref_seq.h"

// Complement base array
const unsigned char comp_base[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32, '!', '"', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/',
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
    '@', 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', '[', '\\',']', '^', '_',
    '`', 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', '{', '|', '}', '~', 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

char ref_antd::get_base(int64_t pos, std::string region) {  // 1-based position
  int len;
  char base = 0;
  if (!region.empty() && this->fai != NULL) {
    seq = fai_fetch(this->fai, region.c_str(), &len);
  }
  if (seq) base = *(seq + (pos - 1));
  free(seq);
  return base;
}

void ref_antd::reverse_complement_codon(char* codon) {
    char temp = comp_base[(unsigned char)codon[2]];
    codon[2] = comp_base[(unsigned char)codon[0]];
    codon[0] = temp;
    codon[1] = comp_base[(unsigned char)codon[1]];
}

char *ref_antd::get_codon(int64_t pos, std::string region,
                          gff3_feature feature) {
  int len;
  seq = fai_fetch(this->fai, region.c_str(), &len);
  int64_t edit_pos = feature.get_edit_position(), codon_start_pos;
  std::string edit_sequence = feature.get_edit_sequence();
  int64_t edit_sequence_size = edit_sequence.size();
  char *codon = new char[3];
  int64_t edit_offset = 0;
  if (pos > edit_pos + edit_sequence_size && edit_pos != -1) {
    edit_offset =
        (pos - edit_pos) > edit_sequence_size
            ? edit_sequence_size
            : (pos - edit_pos);  // Account for edits in position of insertion
  }
  // codon_start_pos is w.r.t the reference sequence
  if (feature.get_strand() == '-'){
    codon_start_pos =
        (feature.get_end() - 1) - feature.get_phase() -
        (((feature.get_end() - feature.get_phase() - (pos + edit_offset))) /
         3) *
            3;
    codon_start_pos -= 2; // This is to get to codon start from the 3' end and then take reverse complement
  } else {
    codon_start_pos =
        (feature.get_start() - 1) + feature.get_phase() +
        (((pos + edit_offset - (feature.get_start() + feature.get_phase()))) /
         3) *
            3;
  }
  for (int i = 0; i < 3; i++) {
    if (codon_start_pos + i < (int32_t)feature.get_start() - 1 ||
        codon_start_pos + i >
            (int32_t)feature.get_end() -
                1) {  // If before or after CDS region return 'N'.
      codon[i] = 'N';
    } else if (codon_start_pos + i < edit_pos - 1 ||
               edit_pos == -1) {  // If before edit or with no edit
      codon[i] = *(seq + codon_start_pos + i);
    } else if (codon_start_pos + i >= edit_pos - 1 &&
               codon_start_pos + i <= edit_pos - 1 + edit_sequence_size -
                                          1) {  // size() - 1 since edit_pos
                                                // include one base already
      codon[i] = edit_sequence[codon_start_pos + i - (edit_pos - 1)];
    } else if (codon_start_pos + i > edit_pos - 1 + edit_sequence_size - 1) {
      edit_offset = (codon_start_pos + i) - (edit_pos - 1) > edit_sequence_size
                        ? edit_sequence_size
                        : (codon_start_pos + i) - (edit_pos - 1);
      codon[i] = *(seq + codon_start_pos + i - edit_offset);
    }
  }

  if (feature.get_strand() == '-') {
      reverse_complement_codon(codon);
  }

  free(seq);
  return codon;
}

char *ref_antd::get_codon(int64_t pos, std::string region, gff3_feature feature,
                          char alt) {
  int len;
  seq = fai_fetch(this->fai, region.c_str(), &len);
  int64_t edit_pos = feature.get_edit_position(), codon_start_pos;
  std::string edit_sequence = feature.get_edit_sequence();
  int64_t edit_sequence_size = edit_sequence.size();
  char *codon = new char[3];
  int i;
  int64_t edit_offset = 0, alt_pos = pos;
  if (pos > edit_pos + edit_sequence_size && edit_pos != -1) {
    edit_offset =
        (pos - edit_pos) > edit_sequence_size
            ? edit_sequence_size
            : (pos - edit_pos);  // Account for edits in position of insertion
  }
  //  TODO: Remove code duplication with function above
  // codon_start_pos is w.r.t the reference sequence
  if (feature.get_strand() == '-'){
    codon_start_pos =
        (feature.get_end() - 1) - feature.get_phase() -
        (((feature.get_end() - feature.get_phase() - (pos + edit_offset))) /
         3) *
            3;
    codon_start_pos -= 2; // This is to get to codon start from the 3' end and then take reverse complement
  } else {
    codon_start_pos =
        (feature.get_start() - 1) + feature.get_phase() +
        (((pos + edit_offset - (feature.get_start() + feature.get_phase()))) /
         3) *
            3;
  }
  for (i = 0; i < 3; ++i) {
    if (codon_start_pos + i < edit_pos - 1 ||
        edit_pos == -1) {  // If before edit or with no edit
      codon[i] = *(seq + codon_start_pos + i);
    } else if (codon_start_pos + i >= edit_pos - 1 &&
               codon_start_pos + i <= edit_pos - 1 + edit_sequence_size -
                                          1) {  // size() - 1 since edit_pos
                                                // include one base already
      codon[i] = edit_sequence[codon_start_pos + i - (edit_pos - 1)];
    } else if (codon_start_pos + i > edit_pos - 1 + edit_sequence_size - 1) {
      edit_offset = (codon_start_pos + i) - (edit_pos - 1) > edit_sequence_size
                        ? edit_sequence_size
                        : (codon_start_pos + i) - (edit_pos - 1);
      codon[i] = *(seq + codon_start_pos + i - edit_offset);
    }
  }
  // Recompute alt position
  edit_offset = 0;
  if (alt_pos > edit_pos + edit_sequence_size && edit_pos != -1) {
    edit_offset =
        (pos - edit_pos) > edit_sequence_size
            ? edit_sequence_size
            : (pos - edit_pos);  // Account for edits in position of insertion
  }
  alt_pos += edit_offset;
  codon[alt_pos - 1 - codon_start_pos] = alt;

  if (feature.get_strand() == '-') {
      reverse_complement_codon(codon);
  }

  free(seq);
  return codon;
}

int ref_antd::add_gff(std::string path) {
  // Read GFF file
  if (!path.empty()) gff.read_file(path);
  return 0;
}

int ref_antd::add_seq(std::string path) {
  this->fai = NULL;
  // Read reference file
  if (!path.empty()) this->fai = fai_load(path.c_str());
  if (!this->fai && !path.empty()) {
    std::cout << "Reference file does not exist at " << path << std::endl;
    return -1;
  }
  return 0;
}

ref_antd::ref_antd(std::string ref_path) {
  this->seq = NULL;
  this->add_seq(ref_path);
}

ref_antd::ref_antd(std::string ref_path, std::string gff_path) {
  this->seq = NULL;
  this->add_seq(ref_path);
  this->add_gff(gff_path);
}

ref_antd::~ref_antd() {
  if (this->fai) fai_destroy(this->fai);
}

// used to add codon info to variants output
int ref_antd::codon_aa_stream(std::string region,
                              std::ostringstream &line_stream,
                              std::ofstream &fout, int64_t pos, char alt) {
  std::vector<gff3_feature> features = gff.query_features(pos, "CDS");
  if (features.size() == 0) {  // No matching CDS
    fout << line_stream.str() << "NA\tNA\tNA\tNA\tNA\tNA" << std::endl;
    return 0;
  }
  std::vector<gff3_feature>::iterator it;
  char *ref_codon, *alt_codon;
  for (it = features.begin(); it != features.end(); it++) {
    fout << line_stream.str();
    // add in gene level info, control for case it's not present
    std::string gene = it->get_attribute("gene");
    if (gene.empty()) {
      fout << it->get_attribute("ID") << "\t";
    } else {
      fout << gene + ":" + it->get_attribute("ID") << "\t";
    }
    ref_codon = this->get_codon(pos, region, *it);
    fout << ref_codon[0] << ref_codon[1] << ref_codon[2] << "\t";
    fout << codon2aa(ref_codon[0], ref_codon[1], ref_codon[2]) << "\t";
    alt_codon = this->get_codon(pos, region, *it, alt);
    fout << alt_codon[0] << alt_codon[1] << alt_codon[2] << "\t";
    fout << codon2aa(alt_codon[0], alt_codon[1], alt_codon[2]) << "\t";

    // adding amino acid position
    // factor in translation direction
    char strand = it->get_strand();
    int64_t aa_pos;
    if (strand == '-') {
      int64_t end = it->get_end();
      aa_pos = ((end - pos) / 3) + 1;
    } else { // when strand is equal to '+', '?', or others
      int64_t start = it->get_start();
      aa_pos = ((pos - start) / 3) + 1;
    }
    fout << aa_pos << std::endl;

    delete[] ref_codon;
    delete[] alt_codon;
  }
  return 0;
}

std::vector<gff3_feature> ref_antd::get_gff_features() {
  return gff.get_features();
}
