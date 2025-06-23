#include "suffix_tree.h"

#include "algorithm"
#include "alignment.h"
#include "fstream"
#include "iostream"
#include "vector"

suffix_node::suffix_node(int b, int *e, suffix_node *p, suffix_node *l) {
  this->children = new suffix_node *[MAX_CHAR];
  for (unsigned int i = 0; i < MAX_CHAR; i++) this->children[i] = 0;
  this->begin = b;
  this->end = e;
  this->nchildren = 0;
  this->parent = p;
  this->suffix_link = l;
}

bool suffix_node::is_leaf_node() { return (nchildren == 0); }

bool suffix_node::contains_child(int ext) {
  return !(this->children[ext] == 0);
}

int suffix_node::get_length() {
  if (this->end == 0) return 0;  // For root
  return *(this->end) - this->begin + 1;
}

int suffix_node::get_depth() {
  int ctr = 0;
  suffix_node *n = this;
  while (n->begin != -1) {
    n = n->parent;
    ctr++;
  }
  return ctr;
}

int get_hamming_distance(std::string s1, std::string s2, int k) {
  int n = 0;
  for (unsigned int i = 0; i < s1.length(); i++) {
    if (s1[i] != s2[i]) {
      n++;
    }
    if (n > k) return -1;
  }
  return n;
}

std::string suffix_node::get_longest_common_substring(std::string s1,
                                                      std::string s2) {
  std::string str = s1 + "#" + s2 + "@";
  std::string p = "", tmp, longest_trvsl = "";
  unsigned int max = 0, i;
  suffix_node *node = this;
  if (node->begin == -1)
    p = "";
  else {
    tmp = node->get_path(str);
    if (s1.find(p + tmp) != std::string::npos &&
        s2.find(p + tmp) != std::string::npos)
      p += tmp;
  }
  tmp = "";
  for (i = 0; i < MAX_CHAR; i++) {
    if (node->children[i] == 0) continue;
    longest_trvsl = node->children[i]->get_longest_common_substring(s1, s2);
    if (longest_trvsl.length() > max) {
      max = longest_trvsl.length();
      tmp = longest_trvsl;
    }
  }
  p += tmp;
  return p;
}

std::string suffix_node::get_path(std::string s) {
  if (this->begin == -1) return "R";
  return s.substr(this->begin, *(this->end) - this->begin + 1);
}

void suffix_node::extend_path(int *e) { this->end = e; }

suffix_node *suffix_node::add_child(int ext, int b, int *e, suffix_node *l) {
  suffix_node *n = new suffix_node(b, e, this, l);
  this->children[ext] = n;
  this->nchildren++;
  return n;
}

void suffix_node::add_child(suffix_node *c, int ext) {
  this->children[ext] = c;
  c->parent = this;
  this->nchildren++;
}

suffix_node *suffix_node::get_child(int ext) { return this->children[ext]; }

bool suffix_node::contains_depth(int depth) {
  return this->get_depth() <= depth &&
         depth <= (this->get_length() + this->get_depth());
}

void suffix_node::print(std::string s) {
  for (int i = 0; i < this->get_depth(); i++) {
    std::cout << " ";
  }
  std::string t =
      (this->begin == -1)
          ? "R"
          : " " + s.substr(this->begin, *(this->end) - this->begin + 1);
  std::cout << t;
  if (this->suffix_link != 0)
    std::cout << " --- ( " << this->suffix_link->parent->get_path(s) << " "
              << this->suffix_link->get_path(s) << ")";
  if (this->begin != -1) std::cout << " - " << this->parent->get_path(s);
  std::cout << std::endl;
  for (unsigned int i = 0; i < alphabet.length(); i++) {
    if (this->children[i] != 0) this->children[i]->print(s);
  }
}

bool suffix_node::walk_next(int &beg, int &suffix_length) {
  suffix_node *node = this;
  if (suffix_length >= node->get_length()) {
    beg += node->get_length();
    suffix_length -= node->get_length();
    return true;
  }
  return false;
}

suffix_node *build_suffix_tree(std::string s) {
  suffix_node *root = new suffix_node(-1, 0, 0, 0);
  root->suffix_link = root;
  root->parent = root;
  suffix_node *cur_node = root, *new_node = 0, *new_cur_node = 0;
  int str_ind[MAX_SIZE], suffix_length = 0, *leaf_end = new int,
                         suffix_count = 0, beg = -1, *end;
  unsigned int i = 0, n_str_ind = 0;
  for (i = 0; i < s.length(); i++) {
    str_ind[i] = alphabet.find(s[i]);
    n_str_ind++;
  }
  *leaf_end = 0;
  // root->add_child(str_ind[0], 0, leaf_end, root);
  for (i = 0; i < n_str_ind; i++) {
    *leaf_end = i;  // Handle Extension Rule 1
    suffix_count++;
    new_node = 0;
    while (suffix_count > 0) {
      if (suffix_length == 0) beg = i;
      if (cur_node->children[str_ind[beg]] == 0) {            // Rule 2
        cur_node->add_child(str_ind[beg], beg, leaf_end, 0);  // Trick 3
        if (new_node != 0) {
          new_node->suffix_link = cur_node;
          new_node = 0;  // No internal node created in Rule 2 here.
        }
      } else {
        new_cur_node = cur_node->children[str_ind[beg]];
        if (new_cur_node->walk_next(beg, suffix_length)) {
          cur_node = new_cur_node;
          continue;
        }
        if (str_ind[new_cur_node->begin + suffix_length] ==
            str_ind[i]) {  // Rule 3
          if (new_node != 0 && cur_node->end != 0) {
            new_node->suffix_link = cur_node;
            new_node = 0;  // No internal node created in Rule 2 here.
          }
          suffix_length++;
          break;
        }
        // Rule 2 - New internal node created
        suffix_node *new_int_node;
        end = new int;
        *end = new_cur_node->begin + suffix_length - 1;
        // New internal node
        new_int_node =
            cur_node->add_child(str_ind[beg], new_cur_node->begin, end, 0);
        new_int_node->add_child(str_ind[i], i, leaf_end, 0);
        new_cur_node->begin += suffix_length;
        new_int_node->add_child(new_cur_node, str_ind[new_cur_node->begin]);
        if (new_node != 0) new_node->suffix_link = new_int_node;
        new_node = new_int_node;
      }
      suffix_count--;
      if (cur_node->end == 0 && suffix_length > 0) {
        suffix_length--;
        beg = i - suffix_count + 1;
      } else if (cur_node->end != 0) {
        if (cur_node->suffix_link != 0) {
          cur_node = cur_node->suffix_link;
        } else {
          cur_node = cur_node->parent->suffix_link;
        }
      } else if (cur_node->get_length() == 1) {
        cur_node = root;
      }
    }
  }
  return root;
}

std::string get_reverse_complement(std::string rev_read) {
  char t = 'N';
  for (unsigned int i = 0; i < rev_read.length(); i++) {
    switch (rev_read[i]) {
      case 'A':
        t = 'T';
        break;
      case 'G':
        t = 'C';
        break;
      case 'C':
        t = 'G';
        break;
      case 'T':
        t = 'A';
        break;
    }
    rev_read[i] = t;
  }
  std::reverse(rev_read.begin(), rev_read.end());
  return rev_read;
}

std::vector<std::string> read_adapters_from_fasta(std::string p) {
  std::vector<std::string> adp;
  std::ifstream fin(p.c_str());
  std::string line;
  while (std::getline(fin, line)) {
    if (line[0] == '>') continue;
    adp.push_back(line);
  }
  return adp;
}

std::string extend_till_mismatches(std::string cmn, std::string l,
                                   std::string adp) {
  if (cmn.length() < MIN_LENGTH) {
    return cmn;
  }
  int pos = adp.find(cmn);
  int r_pos = l.find(cmn);
  if (pos == 0)  // Beginning of adapter match
    return cmn;
  int k = 0;
  pos--;
  r_pos--;
  while (r_pos >= 0 && pos >= 0 && k <= MAX_MISMATCHES) {
    cmn = l[r_pos] + cmn;
    if (l[r_pos] != adp[pos]) k++;
    r_pos--;
    pos--;
  }
  return cmn;
}

int trim_adapter(std::string f1, std::string f2, std::string adp_path,
                 std::string p) {
  std::string l1, l2, l1_rc, l2_rc, s, cmn, trmd_read;
  int *t = new int[2], beg;
  unsigned int pos;
  std::ifstream ff1(f1.c_str()), ff2, pf(p.c_str());
  std::ofstream out((p + ".trimmed.fastq").c_str());
  std::vector<std::string>::iterator it;
  suffix_node *root = (suffix_node *)malloc(sizeof(suffix_node));
  std::vector<std::string> adp = read_adapters_from_fasta(adp_path);
  int i = -1, n = 0, adp_pos;
  if (!f2.empty()) {
    ff2.open(f2.c_str());
    while (std::getline(ff1, l1) && std::getline(ff2, l2)) {
      i++;
      if ((i - 1) % 4 != 0) continue;
    }
  } else {
    while (std::getline(ff1, l1)) {
      i++;
      if ((i - 1) % 4 != 0) {
        out << l1 << std::endl;
        continue;
      }
      beg = 0;
      pos = l1.length();
      trmd_read = l1;
      for (it = adp.begin(); it != adp.end(); ++it) {
        s = l1 + "#" + *it + "@";
        root = build_suffix_tree(s);
        cmn = root->get_longest_common_substring(l1, *it);
        cmn = extend_till_mismatches(cmn, l1, *it);
        if (cmn.empty()) {
          l1_rc = get_reverse_complement(l1);
          s = l1_rc + "#" + *it + "@";
          root = build_suffix_tree(s);
          cmn = root->get_longest_common_substring(l1_rc, *it);
          cmn = extend_till_mismatches(cmn, l1_rc, *it);
          pos = l1_rc.find(cmn);
          beg = l1_rc.length() - pos;
          adp_pos = (*it).find(cmn);
        } else {
          beg = 0;
          pos = l1.find(cmn);
          adp_pos = (*it).find(cmn);
        }
        if (cmn.length() >= MIN_LENGTH && adp_pos != 0) {
          beg = 0;
          pos = l1.length();
          t = align_seqs(l1, *it);
          if (*t != -1) {
            beg = 0;
            pos = *(t + 1);
          } else {
            std::cout << "Reverse!" << std::endl;
            l1_rc = get_reverse_complement(l1);
            t = align_seqs(l1_rc, *it);
            if (*t != -1) {
              pos = l1.length() - *(t + 1);
              beg = *(t + 1);
            }
          }
        }
        std::cout << beg << " " << pos << std::endl;
        if (beg != 0 || pos != l1.length()) {
          trmd_read = l1.substr(beg, pos);
          n++;
          break;
        }
      }
      out << trmd_read << std::endl;
    }
  }
  std::cout << "Number of Adapter Trimmed: :" << n << std::endl;
  delete root;
  delete[] t;
  ff1.close();
  ff2.close();
  pf.close();
  return 0;
}
