#include "parse_gff.h"

gff3_feature::gff3_feature(std::string line) {
  int ctr = 0;
  std::stringstream line_stream;
  std::string cell;
  line_stream << line;
  while (std::getline(line_stream, cell, '\t')) {
    switch (ctr) {
      case 0:  // Name
        this->seqid = cell;
        break;
      case 1:
        this->source = cell;
        break;
      case 2:
        this->type = cell;
        break;
      case 3:
        this->start = atoi(cell.c_str());
        break;
      case 4:
        this->end = atoi(cell.c_str());
        break;
      case 5:
        this->score = atof(cell.c_str());
        break;
      case 6:
        this->strand = cell[0];
        break;
      case 7:
        this->phase = atoi(cell.c_str());
        break;
      case 8:
        this->set_attributes(cell);
        break;
    }
    ctr++;
  }
  if (ctr < 9) std::cout << "GFF file is not in GFF3 file format!" << std::endl;
  line_stream.clear();
}

int gff3_feature::print() {
  std::cout << seqid << "\t" << source << "\t" << type << "\t" << start << "\t"
            << end << "\t" << score << "\t" << strand << "\t" << phase << "\t";
  std::map<std::string, std::string>::iterator it;
  for (it = attributes.begin(); it != attributes.end(); it++) {
    std::cout << it->first << ": " << it->second << "; ";
  }
  std::cout << std::endl;
  return 0;
}

std::string gff3_feature::get_attribute(std::string key) {
  std::string val;
  if (attributes.find(key) != attributes.end()) {
    val = attributes[key];
  }
  return val;
}

int gff3_feature::set_attributes(std::string attr) {
  std::string key, val;
  std::regex exp("[^;]+");
  std::regex_iterator<std::string::iterator> it(attr.begin(), attr.end(), exp);
  std::regex_iterator<std::string::iterator> rend;
  std::string delimiter = "=";
  while (it != rend) {
    key = it->str().substr(0, it->str().find(delimiter));
    val = it->str().substr(it->str().find(delimiter) + 1, it->str().length());
    if (!key.empty() && !val.empty()) {
      this->attributes[key] = val;
    }
    ++it;
  }
  return 0;
}

uint64_t gff3_feature::get_start() { return start; }

uint64_t gff3_feature::get_end() { return end; }

int gff3_feature::get_phase() { return phase; }

std::string gff3_feature::get_type() { return type; }

int64_t gff3_feature::get_edit_position() {
  int64_t edit_pos = -1;
  std::map<std::string, std::string>::iterator it;
  for (it = attributes.begin(); it != attributes.end(); it++) {
    if (it->first.compare(EDIT_POSITION) == 0) {
      edit_pos = stoi(it->second);
      break;
    }
  }
  return edit_pos;
}

std::string gff3_feature::get_edit_sequence() {
  std::string edit_seq = "";
  std::map<std::string, std::string>::iterator it;
  for (it = attributes.begin(); it != attributes.end(); it++) {
    if (it->first.compare(EDIT_SEQUENCE) == 0) {
      edit_seq = it->second;
      break;
    }
  }
  return edit_seq;
}

int gff3::print() {
  std::vector<gff3_feature>::iterator it;
  for (it = features.begin(); it != features.end(); it++) {
    it->print();
  }
  return 0;
}

std::vector<gff3_feature> gff3::get_features() { return features; }

gff3::gff3() { this->is_empty = true; }

gff3::gff3(std::string path) {
  this->is_empty = true;
  this->read_file(path);
}

int gff3::read_file(std::string path) {
  std::ifstream fin = std::ifstream(path);
  if (!fin) {
    std::cout << "GFF file does not exist at " << path << std::endl;
    return -1;
  }
  std::string line;
  while (std::getline(fin, line)) {
    if (line[0] == '#' || line.empty())  // Avoid comments in GFF file
      continue;
    features.push_back(gff3_feature(line));
  }
  if (this->get_count() > 0) this->is_empty = false;
  return 0;
}

std::vector<gff3_feature> gff3::query_features(uint64_t pos, std::string type) {
  std::vector<gff3_feature>::iterator it;
  std::vector<gff3_feature> res;
  for (it = features.begin(); it != features.end(); it++) {
    if (it->get_type() != type) continue;
    if (pos >= it->get_start() && pos <= it->get_end()) {
      res.push_back(*it);
    }
  }
  return res;
}

int gff3::get_count() { return features.size(); }

bool gff3::empty() { return is_empty; }
