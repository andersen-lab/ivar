#include <iostream>
#include <vector>

#include "../src/parse_gff.h"

int main() {
  int success = 0;

  gff3 *gff = new gff3();
  success += gff->read_file("./data/test.gff");
  success += gff->empty();
  return (success == 0);
}
