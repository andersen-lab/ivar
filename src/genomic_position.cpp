#include <string>
#include <vector>

void position::update_alleles(std::string nt, uint32_t count, uint32_t qual){
  //update overall positions depth
  if(nt.find("+") == std::string::npos){
    depth += count;
  }
  //check if in allele vector
  int exists = check_allele_exists(nt, alleles);
  //allele does not exist
  if (exists == -1){
    allele tmp;
    tmp.mean_qual = qual;
    tmp.depth = count;
    tmp.nuc = nt;
    alleles.push_back(tmp);
  } else {
    alleles[exists].mean_qual += qual;
    alleles[exists].depth += count;
  }
}

int check_position_exists(uint32_t p, std::vector<position> positions) {
  for (uint32_t i=0; i < positions.size(); i++) {
    if (p == positions[i].pos) {
      return((int)i);
    }
  }
  return(-1);
}
