#include <iostream>

#include "../src/allele_functions.h"

int main() {
  int num_tests = 7;
  std::string b = "AAAATTTG+3ATGT-3ATG";
  int size = 9, success = 0, i;
  char ref = 'A';
  uint8_t _q[] = {30, 30, 20, 23, 20, 15, 20, 40, 20};
  std::string q = "";
  for (int i = 0; i < size; ++i) {
    _q[i] += 33;
    q += _q[i];
  }
  std::cout << q << std::endl;
  std::vector<allele> ad = update_allele_depth(ref, b, q, 10);
  i = find_ref_in_allele(ad, 'T');
  success += (ad.at(i).mean_qual == 18) ? 1 : 0;
  i = find_ref_in_allele(ad, 'A');
  success += (ad.at(i).mean_qual == 25) ? 1 : 0;

  ad.clear();
  ad = update_allele_depth(ref, b, q, 25);
  success += (ad.size() == 4) ? 1 : 0;
  i = find_ref_in_allele(ad, 'A');
  success += (ad.at(i).mean_qual == 30) ? 1 : 0;
  i = find_ref_in_allele(ad, 'G');
  success += (ad.at(i).mean_qual == 40) ? 1 : 0;

  ref = 'T';
  b = "TtTtTTTtTttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "ttttttttgttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttattttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttatttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "ttagtttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt"
      "ttttttttttttttttttttttttttt*";
  q = "GImGmo.IGGIGGGGGCGEGGI9GGGGGGIGIIGGAIGIGIGGIGGGFDGDGGGFGGEGGGEEGGGGGIG."
      "IIIGIIIGIIIGIAIAIAGGIIAG.,GGGGG.IIIGIIGGIGIIGGGIGGGIFG="
      "F9GGGGGDGGGGGGEGFGGDG;FDGG@GGGC?GGGGG@GGGGGGFFFGGGFGG;G<,"
      "GGGFGGGCGGGG9GGAF9GGGGGGDGGAIGIIIGIGIGGIIIIGAIIGIIFGFFGDGGGCGGCGFGGGG9GG"
      "GGFGDDFGFGGGGFGFGGGGGFGGGGGFGFGFAGGFFGGGG=GCGEGGEGEGG9C?"
      "GGFEGGGGGGGGGGGG9?CGGGGEC;GGGGFEEGGGGFFFGGGGGEGFGEF=GG9GDF<"
      "FGGGGFGGGGFGGFAA4D?G9GDGGG@FGGGGGGGGF,GGGGGGGGGGG,GFCB;"
      "CGGGGDGGGGGGGGGGFGAFGGGE@F,G<GFGGGGGDGGGGGGG@GG8GFDGGGGGGGGCGFAGGGGGGGG="
      "GGGGGGDFGGFGGGGGFGGGGGGG@GEGEEGGGGGEGGGFFGGGGGFGGGGG;"
      "GGDGGGGGGGGGFGGGGFG,GGGGGGGGFGG9E@GGGFGGGDGGFFGG4G,"
      "GCGFG9GFFGGGGGFGGGGGGGCGGGG6GGFGGEGGF<GFGGGGGGGA9GGGG,E,FGGG,GFGG<BG<"
      "GGEGGGGG,"
      "BG9EFGGAGGFFFGFGGFGGFFGGEFEGGAGGGGGGGGGG9G9GDGGGGGGGGGGGGGGFGGGGGFCGGGGG"
      "GFGGGGGD@FGGFFGGEGG@GCGGGGFDGGGFG9G9GGGGAGCGGGG<"
      "GGFGGGCGGGGGGGGGGGGGFGGFGGFGGGGGGGFGGFEG,GFGFGGG?FCGG=GCF:"
      "FGCEGCGFGGDGGGGGDGGGGGGGGGGGGGGGGGFGDGGGAGGGG@GGG=GFGGGCGECCAFFE:"
      "GGGGGGGFFG<"
      "FGGGCGGGGGGGGGGGCGGGGGFFFGFG9GGGGGGGGGG3GGGGGGGGFGFG5GGGGGGGGGGGGGGGF7FF"
      "GFFGF8GGGGGGGFGFG;,GGFG,FBGFGGGGFGGGGEGEGGGGG9GGEFGGFCGGFFFGAGFGGGGGGF<"
      "GFEFGFGGGCGGFEFG<EFGGCGGFGAGGGEACFFGGGGFGFGGGE9GGGGCG,G9GG9GGG@GCGGG<"
      "GGGGGFEGGGF,GGFD@GFCFGG;G?FGGEGGGGGGGFG=GFGGGGG9DGFGGGGEGGFCGFGGGGEGDG<"
      "GGFGFGGGGGGGCFGEGDGG,G9GAG,GGEFEGGG,C>"
      "GGGGGGGGGGGGGGGGGCGG5GCGDGGGG9GGGGGGGGGGGGGGGGFGGG=GGGGE>CCGGEGGGAFGFD?"
      "GCG9FGAGG@GGGGGAG,ADGEGGGGEGGGGGGFGGGGCGGGGGFEEDGGGGGDFG9GGEGGFFC@"
      "DGGGGGGGFGGGFFGGGGGEGFGGGGGGGGGGGGDGGGGGGFGG<GGGFAGGGGGGGGGGEGG?"
      "GGGGGG9GGCGAGG,GFGGGE,FGG9FGGGGGGGCGGCD99GGF,GGGGGGGGFFGGG;"
      "AGG9GGGGGG9GDGG9GFGDGGDGGGCG9GG?GGAGAAGGGG="
      "9GAGCGGGGGGGGGGGGGGGGGGGGGGCGFGFFGGGGGGGEGGGGG,,GGGGGFCGDDF;"
      "GGGCGGFGEGGDGGGG,"
      "GGFGGGGFGGGGGGGGCGGGGGGGGGGFGGGFGGGGGDCCGGGGGGFGGGFGGGGGGGGGGGG?FAGG;"
      "GGGCGGGGFGG,GDGFFGGEGGFGCGGFGGFC>GFGGFAG?;,DGGGFGGGFGCGAGEG@GGGFGGG,"
      "GG9GDFFFFAGFGFCGGGGF9?GCGGGGGGGFGFGGGGFG@FGGFGGGGFFGGGGGFGGFGGFG:;GGGG@"
      "FGGFGGEEFDGGGGFGEGGD9GFFEDGGAGGGGGFCGGCGGF;"
      "GF9GGGFFGAAG9GEGGGGGGBFGGGGEFFGCGGCGFE?GGG9;5CGFGFDFG9GAGGGGGGGGGG<"
      "GGGGGGFGEGGGGF?;GGGGGAGGGFGGEGG;GGG;GGFG,GGGGGGGGCGFGGGGG;"
      "GGGGAGGGFGGGGGFDG,CGDFGFGG<GGAGGFGE9FG,;FFGGFFGGGEGGFGGGEFGCGGFGGG,G<,"
      "GFGGFGGFGGGGGGDGCEGEGGGGGGGFFGGGF4GGFAGGEG9GGGGGBGF9GEAGEG9GG9GGGGFGGGGG"
      "GGGFG,GGAFGAFGGGGG<GGG;GGAGGGBGGGFGFGGFGGGGGG@GGGEFG9FFEG,"
      "GGDGGGGGGGGGGGDGGGFGGGFGGGGGGCGGGGCGFGGGGG9GGFFCGCEGGG,G;"
      "GFGAGFCGFFGGGGGGGEGGGFGGFFFGGGGGGGGGGGFGGGFCFGG,GGCGGGGGGGFAGGFAGGGGE?"
      "GF,AGGGGBDGGGG@CGGCG9EGGBGGEGEGGGCGECAGGGFGGFG,DGGC,9E9GGGC,FFG,FCGDGG,"
      "GGCFG;GGGG,GGGGCGG@,G9GGGGGGGGGFFAGGG,"
      "AFGGGGGGCFGGGFFFGFGGCGGGCFGGAEFEGGGGGGGDGGGGGGGF9GGGDGGCGGGGGGGGGGGGGGFG"
      "DGGFG?GDCE5EGGGGGGGF,FGEGEG,GGEGGGGGGGGEG4FGGGGGGGGG?GF,"
      "GGG9FFGGGGGGCCGGGG,GGGCCFE@GGFGDGFGEGGF,FGG,"
      "DEGGGGGGGGGGGGGGGDGGGFAGGFGFGGGGEGGGGFFGAGGGFFFFGGGGGGGG="
      "FGGFFGGGGGGGGGGGCGGFGGGGGGFGGGG<GGGGCCGGGGEFG4GGGGGGFGCGGCFCF,GDDGG?"
      "FGGEEGGEFGGGGD<GGEGGGGGGCGFGFECGFGFGGGGEGGGGGGGGGE9GGGFGGFGFC@GEF@G,,"
      "FGGGGGEGG=GGFCEFGGFGGGFGEGGGGGGGGGGGEFGCGFGGFGGGGGGG;G,GGGG?"
      "DAGGFFFG4GGEFGEFGCGGFGGFG,GGGGGGGGGGGAG9GGGGGGGGGGGFGGFC;"
      "GGFGFGGGGGG4GFGEGGGFEGGGGG:G;GG9GCGAG,FGGGFEGFGGGG,"
      "GGGGGGDGGGGGGGGDF9CCG;GGFGFGFGDDGGGFGCGCGFGDGDFGGFFIGGIGIAGGIGIGI."
      "IIIIII";

  std::cout << std::endl << std::endl << std::endl;
  ad = update_allele_depth(ref, b, q, 20);
  i = find_ref_in_allele(ad, 'T');
  success += (ad.at(i).mean_qual == 36) ? 1 : 0;
  ref = 'A';
  b = "...............................,,,,,,,,,,,,,,,,t,,,,,,,,,..-1T......-"
      "2TG....................................................................."
      "........................................................................"
      "........................................................................"
      "...................................................T...................."
      "............................................-1T........................."
      "........................................................................"
      "......................................................-1T..............."
      "........................................................................"
      ".......,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,+2tg,,,,,,,,,,,,,,,,,"
      ",,,,,,,,,,,,,,,,,,,,,,,+2tg,,,,,,,,,,,,,,,,,,,,,,+2tg,,,-1t,,,,,,,,,,,,,"
      ",,,,,,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,"
      "t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,-16tggaaaaaggagtggt,,,,"
      ",,,,,,,,,,,,,,,-1t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,................"
      "........................................................................"
      "........................................................................"
      "......................................................T................."
      ".................T......................................................"
      "........................................................................"
      "........................................................................"
      "........G.......G......................................................."
      "..................................................................,,,,,,"
      ",t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
      ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,-1t,,,,,,,,,,,,,,,,,,,,,,,,,,,"
      ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,"
      ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,"
      ",,,,,t,,.";
  q = "../2259/8028.4931;44//029.45090../213..03303..2/3//00212/.011/4./27/"
      "9.13005.18211::.2/81//25351/23351717/0222/:2271//.3..2647.9123//0/"
      "3.56:9/75/:4/2371/01:;6721289/.138.02/.640231/710.37:90/16/91101482.6.1/"
      ".00/13.9;/.010.232813310/29.5120130610/1.970172.60/24/11/204:60..7/120//"
      "07/..91.//432.1.3/3.:82.1.525676247/;1.;0:40.2184.23364.0/6:.053410320/"
      "91413/01:765.1.8//1/00..2/::515/.162//42932/"
      "0011240346394013:20442017;.121.3601419.31.2199900106.332/2/"
      "3.322.0114.:939/899/208/1;0090/28/.:3.1:3728690/2:19.30.030../1:/:1./"
      "397.9;./00/0402.299825/1:/.6.71/..0012384./431/6/20//127/..5/7246/"
      "000:5841650.03/38550..90125822.753131.6994085/6141.8...1.//.1.1//4.1.//"
      "...//0.9./.2../00201.2....//6/./0/5/0002/.3/01/.0....2.0/.0.0.00/81/"
      ".010.0.41././///0.301/410/./00.1..60/0.2/.0311.0//...//2.1./2/0./0//./"
      "..0/1.1.../2./7/3.000/...432.0.5.1./000.200/.20.2..131.0/3..0..///04/"
      "00.//3..40..4./1/40001.5./400//...5.0/0.4....007./..705.0.10/2.0...//"
      "11..108357/36125/5.103.7.52110210111706.4.8/1111:./00.925.3/21;01.1/"
      "0.3.90/.32.41741//6/..1760871:.04.20/.29.1.19.310/6709001330/0/.3/95.42/"
      "00.19820//.05/717631./7/50:/.4151.8..01.4:980.4/5/"
      "60.02751241140210940:2/7579/.2.1.9:40/0/21013/2./8/.0///;7:12132..00./"
      "64/0/2:222410/3/40.9..1.491.1/1:9/11010/01/02.0004:/.9632.216//320/.0/"
      "3.0/7.171/.6922030133030.0323.97802511151/..3:/5973.7/.//080.11.230.801/"
      ".:1.2.9133/11.0.2200//.6.24/:202971.12.1/./042/16310/23881.8///"
      "00.003;2;:/91:6601.64.71706107..2/"
      "4241.003.7116.:9009.000890952.24.032.90//1.26102343/31928://"
      "995010.11.55/72.123:6011:117/:0.061.0/1.1//0130/...3/82/002/1./.2/03/0/"
      "./...3./1/.34./0/000..2..01.../11./.10/.120/0.///11.0.//0.10.0121...3./"
      "000./5/33////0//0/0/.0..0/.0//0../3.03201/0/.0//.0/1/.6.5.201/102.0.020/"
      "//3..11/111.11.///.0202/../20../.10//8.0.//.1./.0./.8/..////.//"
      "30...003.1/0/0/12../1140.2.1010./11./3.3/025/../0.0/18";
  ad = update_allele_depth(ref, b, q, 10);
  print_allele_depths(ad);
  i = find_ref_in_allele(ad, 'T');
  success += (ad.at(i).depth == 12) ? 1 : 0;
  return (num_tests == success) ? 0 : -1;
}
