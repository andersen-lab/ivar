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

  if(!g[2].get_previous())
    num_success += 1;
  if(g[2].get_next())
    num_success += 1;
  if(g[2].get_next()->get_start() == 203)
    num_success += 1;
  if(g[3].get_previous()->get_end() == 26)
    num_success += 1;

//  char *codon = new char[3];
//  codon = refantd.get_codon(30, "test", g.at(2));
//  num_success = check_failure('C', 'A', 'T', codon, 30);
//  std::cout << num_success;
  if (num_success == 3) return 0;
  return -1;
}
