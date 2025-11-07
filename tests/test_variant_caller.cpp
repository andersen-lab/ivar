#include "../src/variant_caller.h"

void print_cigar(uint32_t *cigar, int nlength) {
  for (int i = 0; i < nlength; ++i) {
    std::cerr << ((cigar[i]) & BAM_CIGAR_MASK);
    std::cerr << "-" << ((cigar[i]) >> BAM_CIGAR_SHIFT) << " ";
  }
  std::cerr << std::endl;
}

struct expected_read {
  std::vector<uint32_t> positions;
  std::vector<std::string> bases;
  std::vector<uint32_t> qualities;
};

std::vector<expected_read> expected_reads = {
    {{9, 42}, {"-C", "G"}, {20, 52}}, // Read name: test_outside_primer
    {{30, 42}, {"+G", "T"}, {63, 49}}, // Read name: test_20_381_1:0:0_6:0:0_1_150
    {{20, 25, 23}, {"C", "+A", "-CT"}, {20, 63, 20}}, // Read name: test_20_381_1_softclip:0:0_6:0:0_1_150
    {{353, 365}, {"-ACCA", "T"}, {20, 27}} // Read name: test_20_381_1:0:0_6:0:0_1_150
};

int main(){
  uint8_t min_qual = 20;
  std::string ref_path="../../data/db/test_ref.fa";
  variant_caller vc(min_qual, ref_path, "");
  samFile *bam = sam_open("../../data/test.unmapped.sorted.bam", "r");
  sam_hdr_t *header = sam_hdr_read(bam);
  bam1_t *aln = bam_init1();
  int ctr = 0, test_idx = 0;
  std::vector<uint32_t> positions1 = {};
  std::vector<std::string> bases1;
  std::vector<uint32_t> qualities1;

  while (sam_read1(bam, header, aln) >= 0) {
    if (ctr != 0 && ctr != 1 && ctr != 2 && ctr != 7) {
      ctr ++;
      continue;
    }
    std::vector<uint32_t> positions;
    std::vector<std::string> bases;
    std::vector<uint32_t> qualities;

    vc.parse_read(aln, "test", positions, bases, qualities);

    for(int i = 0; i < expected_reads[test_idx].positions.size(); i++){
      int test_position = expected_reads[test_idx].positions[i];
      std::vector<int> test_indices;
      // For insertions, there may be multiple positions matching
      for (size_t j = 0; j < positions.size(); ++j) {
        if (positions[j] == test_position) {
          test_indices.push_back(j);
        }
      }
      if(test_indices.empty()){
        std::cerr << "Could not find position " << test_position << " in read " << bam_get_qname(aln) << "\n";
        return 1;
      }
      bool pass = false;
      for (int index : test_indices) {
        if(bases[index] == expected_reads[test_idx].bases[i] && qualities[index] == expected_reads[test_idx].qualities[i]){
          pass = true;
        }
      }
      if(!pass){
        std::cerr << "Position mismatch at read " << bam_get_qname(aln) << ": at position " << test_position << " expected " << expected_reads[test_idx].bases[i] << "(" << expected_reads[test_idx].qualities[i] << ")" << " got\n";
        for(int index: test_indices){
          std::cerr << bases[index] << "(" << qualities[index] << ")" << "\n";
        }
        return 1;
      }
    }

    ctr++;
    test_idx++;
  }
  return 0;
}