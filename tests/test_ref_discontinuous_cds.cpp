#include "../src/ref_seq.h"

int main() {
  int success = 0;
  ref_antd refantd("../../data/db/test_ref.fa", "../../data/test_discontinuous_cds.gff");

  // + strand
  std::vector<cds_group> fwd_groups = refantd.query_cds_groups(15);
  success += fwd_groups.size() == 1;
  if(success != 1) {
    std::cout << "Error querying forward group at pos 15";
    return -1;
  }

  // pos 10: first codon in seg 1, cds_pos 0, expect AGC
  char *codon = refantd.get_codon(10, "test", fwd_groups[0]);
  success += (codon[0] == 'A' && codon[1] == 'G' && codon[2] == 'C');
  if(success != 2) {
    std::cout << "Error at fwd pos 10 (expected AGC): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 29: cds_pos 19, codon on boundary -> ref 28,29,100 = TCT
  codon = refantd.get_codon(29, "test", fwd_groups[0]);
  success += (codon[0] == 'T' && codon[1] == 'C' && codon[2] == 'T');
  if(success != 3) {
    std::cout << "Error at fwd pos 29 (straddle, expected TCT): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 100: cds_pos 20, same codon TCT (snps on previous CDS)
  codon = refantd.get_codon(100, "test", fwd_groups[0]);
  success += (codon[0] == 'T' && codon[1] == 'C' && codon[2] == 'T');
  if(success != 4) {
    std::cout << "Error at fwd pos 100 (straddle, expected TCT): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 101: cds_pos 21, first codon entirely in seg 2 -> ref 101,102,103 = CTC
  codon = refantd.get_codon(101, "test", fwd_groups[0]);
  success += (codon[0] == 'C' && codon[1] == 'T' && codon[2] == 'C');
  if(success != 5) {
    std::cout << "Error at fwd pos 101 (expected CTC): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 50: outside the joined CDS -> NNN
  codon = refantd.get_codon(50, "test", fwd_groups[0]);
  success += (codon[0] == 'N' && codon[1] == 'N' && codon[2] == 'N');
  if(success != 6) {
    std::cout << "Error at fwd pos 50 (out of CDS, expected NNN): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // '-' strand'
  std::vector<cds_group> rev_groups = refantd.query_cds_groups(319);
  success += rev_groups.size() == 1;
  if(success != 7) {
    std::cout << "Error querying reverse group at pos 319";
    return -1;
  }

  // pos 319: cds_pos 0, top=CAG (5' to 3'), after complement = GTC
  codon = refantd.get_codon(319, "test", rev_groups[0]);
  success += (codon[0] == 'G' && codon[1] == 'T' && codon[2] == 'C');
  if(success != 8) {
    std::cout << "Error at rev pos 319 (expected GTC): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 220: cds_pos 20, top=ACA (ref 301,300,220 in 5' to 3'),
  // after complement = TGT straddles CDSs
  codon = refantd.get_codon(220, "test", rev_groups[0]);
  success += (codon[0] == 'T' && codon[1] == 'G' && codon[2] == 'T');
  if(success != 9) {
    std::cout << "Error at rev pos 220 (straddle, expected TGT): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // + strand within-group overlap: id-overlap-rev = [40,50] U [50,59].
  std::vector<cds_group> fwd_overlap_groups = refantd.query_cds_groups(50);
  success += fwd_overlap_groups.size() == 1;
  if(success != 10) {
    std::cout << "Error querying fwd-overlap group at pos 50";
    return -1;
  }

  // pos 50: boundary base. genomic_to_cds_pos returns first hit (cds 10).
  // codon 3 -> reads cds 9,10,11 = genomic 49,50,50 -> ref CTT.
  codon = refantd.get_codon(50, "test", fwd_overlap_groups[0]);
  success += (codon[0] == 'C' && codon[1] == 'T' && codon[2] == 'T');
  if(success != 11) {
    std::cout << "Error at fwd-overlap pos 50 (boundary, expected CTT): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 52: cds 12, codon 4 entirely in second segment -> 51,52,53 = TGG.
  codon = refantd.get_codon(52, "test", fwd_overlap_groups[0]);
  success += (codon[0] == 'T' && codon[1] == 'G' && codon[2] == 'G');
  if(success != 12) {
    std::cout << "Error at fwd-overlap pos 52 (expected TGG): " << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // - strand within-group overlap: id-overlap-rev = [160,170] U [170,179].
  std::vector<cds_group> rev_overlap_groups = refantd.query_cds_groups(170);
  success += rev_overlap_groups.size() == 1;
  if(success != 13) {
    std::cout << "Error querying rev-overlap group at pos 170";
    return -1;
  }

  // pos 170: boundary base. genomic_to_cds_pos returns first hit (cds 9).
  // codon 3 -> reads cds 9,10,11 = genomic 170,170,169 = ref CCT -> complement = GGA.
  codon = refantd.get_codon(170, "test", rev_overlap_groups[0]);
  success += (codon[0] == 'G' && codon[1] == 'G' && codon[2] == 'A');
  if(success != 14) {
    std::cout << "Error at rev-overlap pos 170 (boundary, expected GGA): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 175: cds 4, codon 1 first sorted segment -> 176,175,174
  // = ref GTC -> complement = CAG.
  codon = refantd.get_codon(175, "test", rev_overlap_groups[0]);
  success += (codon[0] == 'C' && codon[1] == 'A' && codon[2] == 'G');
  if(success != 15) {
    std::cout << "Error at rev-overlap pos 175 (expected CAG): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  // pos 165: cds 15, codon 5 entirely in second sorted segment -> 165,164,163
  // = ref AAG -> complement = TTC.
  codon = refantd.get_codon(165, "test", rev_overlap_groups[0]);
  success += (codon[0] == 'T' && codon[1] == 'T' && codon[2] == 'C');
  if(success != 16) {
    std::cout << "Error at rev-overlap pos 165 (expected TTC): "
              << codon[0] << codon[1] << codon[2];
    delete[] codon;
    return -1;
  }
  delete[] codon;

  return 0;
}
