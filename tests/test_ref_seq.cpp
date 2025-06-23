#include "../src/ref_seq.h"

int check_failure(char n1, char n2, char n3, char *codon, int pos) {
  if (codon[0] != n1 || codon[1] != n2 || codon[2] != n3) {
    std::cout << "Pos: " << pos << " "
              << "Codon: " << n1 << n2 << n3 << " Res: " << codon << std::endl;
    return -1;
  }
  return 0;
}

int main() {
  int num_success = 0;
  ref_antd refantd("../data/db/test_ref.fa", "../data/test.gff");
  std::vector<gff3_feature> g = refantd.get_gff_features();
  char *codon = new char[3];
  codon = refantd.get_codon(30, "test", g.at(2));
  num_success = check_failure('C', 'A', 'T', codon, 30);
  codon = refantd.get_codon(30, "test", g.at(3));
  num_success = check_failure('T', 'C', 'A', codon, 30);
  codon = refantd.get_codon(29, "test", g.at(2));
  num_success = check_failure('C', 'A', 'T', codon, 29);
  // After edit_pos
  codon = refantd.get_codon(107, "test", g.at(4));
  num_success = check_failure('C', 'A', 'A', codon, 107);
  codon = refantd.get_codon(105, "test", g.at(4));
  num_success = check_failure('C', 'A', 'T', codon, 105);
  codon = refantd.get_codon(101, "test", g.at(4));
  num_success = check_failure('T', 'C', 'T', codon, 101);
  codon = refantd.get_codon(100, "test", g.at(4));
  num_success = check_failure('G', 'G', 'A', codon, 100);
  // 2 AA insertion
  codon = refantd.get_codon(100, "test", g.at(5));
  num_success = check_failure('G', 'G', 'T', codon, 100);
  codon = refantd.get_codon(101, "test", g.at(5));
  num_success = check_failure('C', 'A', 'A', codon, 101);
  codon = refantd.get_codon(140, "test", g.at(5));
  num_success = check_failure('C', 'T', 'A', codon, 140);
  std::cout << num_success;
  if (num_success == 0) return 0;
  return -1;
}
