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
  ref_antd refantd("../data/db/test_ref.fa", "../data/test_multiple_segments.gff");
  std::vector<gff3_feature> g = refantd.get_gff_features();
  num_success += g.at(2).get_start() == 1;

//  char *codon = new char[3];
//  codon = refantd.get_codon(30, "test", g.at(2));
//  num_success = check_failure('C', 'A', 'T', codon, 30);
//  std::cout << num_success;
  if (num_success == 1) return 0;
  return -1;
}
