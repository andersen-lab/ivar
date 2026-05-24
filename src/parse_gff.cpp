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
  if (ctr < 9) std::cerr << "GFF file is not in GFF3 file format!" << std::endl;
  line_stream.clear();
}

int gff3_feature::print() {
  std::cerr << seqid << "\t" << source << "\t" << type << "\t" << start << "\t"
            << end << "\t" << score << "\t" << strand << "\t" << phase << "\t";
  std::map<std::string, std::string>::iterator it;
  for (it = attributes.begin(); it != attributes.end(); it++) {
    std::cerr << it->first << ": " << it->second << "; ";
  }
  std::cerr << std::endl;
  return 0;
}

std::string gff3_feature::get_attribute(std::string key) const {
  std::string val;
  std::map<std::string, std::string>::const_iterator it = attributes.find(key);
  if (it != attributes.end()) {
    val = it->second;
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

uint64_t gff3_feature::get_start() const { return start; }

uint64_t gff3_feature::get_end() const { return end; }

int gff3_feature::get_phase() const { return phase; }

std::string gff3_feature::get_type() const { return type; }

char gff3_feature::get_strand() const { return strand; }

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
    std::cerr << "GFF file does not exist at " << path << std::endl;
    return -1;
  }
  std::string line;
  while (std::getline(fin, line)) {
    if (line[0] == '#')  // Avoid comments in GFF file
      continue;
    if(!line.empty())
      features.push_back(gff3_feature(line));
  }
  if(!features.empty()){
    this->is_empty = false;
  } else {
    std::cerr << "GFF file is empty!" << std::endl;
  }

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

cds_group::cds_group() : strand_('+') {}

void cds_group::add_segment(const gff3_feature &f) {
  if (segments_.empty()) {
    strand_ = f.get_strand();
    id_ = f.get_attribute("ID");
    if (id_.empty()) id_ = f.get_attribute("Parent");
    gene_ = f.get_attribute("gene");
  }
  segments_.push_back(f);
}

void cds_group::sort_and_finalize_segments() {
  if (segments_.empty()) return;
  if (strand_ == '-') {
    std::sort(segments_.begin(), segments_.end(),
              [](const gff3_feature &a, const gff3_feature &b) {
                return a.get_start() > b.get_start();
              });
  } else {
    std::sort(segments_.begin(), segments_.end(),
              [](const gff3_feature &a, const gff3_feature &b) {
                return a.get_start() < b.get_start();
              });
  }
  cumulative_len_before_.assign(segments_.size(), 0);
  for (size_t i = 1; i < segments_.size(); ++i) {
    int64_t prev_len =
        segments_[i - 1].get_end() - segments_[i - 1].get_start() + 1;
    cumulative_len_before_[i] = cumulative_len_before_[i - 1] + prev_len;
  }
}

int64_t cds_group::genomic_to_cds_pos(int64_t pos) const {
  for (size_t i = 0; i < segments_.size(); ++i) {
    int64_t start = segments_[i].get_start();
    int64_t end = segments_[i].get_end();
    if (pos < start || pos > end) continue;
    if (strand_ == '-') {
      return cumulative_len_before_[i] + (end - pos);
    }
    return cumulative_len_before_[i] + (pos - start);
  }
  return -1;
}

int64_t cds_group::cds_to_genomic_pos(int64_t cds_pos) const {
  if (cds_pos < 0) return -1;
  for (size_t i = 0; i < segments_.size(); ++i) {
    int64_t seg_len =
        segments_[i].get_end() - segments_[i].get_start() + 1;
    int64_t lo = cumulative_len_before_[i];
    int64_t hi = lo + seg_len;
    if (cds_pos < lo || cds_pos >= hi) continue;
    int64_t offset = cds_pos - lo;
    if (strand_ == '-') {
      return segments_[i].get_end() - offset;
    }
    return segments_[i].get_start() + offset;
  }
  return -1;
}

char cds_group::get_strand() const { return strand_; }

int cds_group::get_phase() const {
  return segments_.empty() ? 0 : segments_[0].get_phase();
}

int64_t cds_group::length() const {
  int64_t total = 0;
  for (size_t i = 0; i < segments_.size(); ++i) {
    total += segments_[i].get_end() - segments_[i].get_start() + 1;
  }
  return total;
}

bool cds_group::contains(int64_t pos) const {
  for (size_t i = 0; i < segments_.size(); ++i) {
    int64_t start = segments_[i].get_start();
    int64_t end = segments_[i].get_end();
    if (pos >= start && pos <= end) return true;
  }
  return false;
}

const std::string &cds_group::get_id() const { return id_; }
const std::string &cds_group::get_gene() const { return gene_; }
const std::vector<gff3_feature> &cds_group::segments() const {
  return segments_;
}
