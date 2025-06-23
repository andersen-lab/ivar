#include "alignment.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

/* Substitution Matrix

   A  T  G  C N
A  1 -1 -1 -1 0
T -1  1 -1 -1 0
G -1 -1  1 -1 0
C -1 -1 -1  1 0
N  0  0  0  0 0
*/

int get_sub_score(char a, char b) {
  int i = 4, j = 4;
  switch (a) {
    case 'A':
      i = 0;
      break;
    case 'T':
      i = 1;
      break;
    case 'G':
      i = 2;
      break;
    case 'C':
      i = 3;
      break;
    case 'N':
      i = 4;
      break;
  }
  switch (b) {
    case 'A':
      j = 0;
      break;
    case 'T':
      j = 1;
      break;
    case 'G':
      j = 2;
      break;
    case 'C':
      j = 3;
      break;
    case 'N':
      j = 4;
      break;
  }
  return substitution[i][j];
}

int get_gap_penalty(int k, char a) {
  if (a == 'N')  // Allow for i7 and i5 indexesto be specified as NNNNNNNN
    return 0;
  return k * gap_open;  // Linear Gap Penalty
}

void print_matrix(int h[][max_adapter_size], int r, int c, std::string read,
                  std::string adap) {
  for (int i = 0; i < r; i++) {
    if (i == 0) {
      std::cout << "    ";
      for (int k = 0; k < c - 1; k++) std::cout << adap[k] << " ";
      std::cout << "\n";
    }
    for (int j = 0; j < c; j++) {
      if (j == 0 && i > 0) {
        if (i - 1 < (int)read.length()) std::cout << read[i - 1] << " ";
      }
      if (i == 0 && j == 0)
        std::cout << " "
                  << " ";
      std::cout << h[i][j] << " ";
    }
    std::cout << "\n";
  }
}

// Return Value 1 - Diag, 2->Right, 3->Down, 0-> 0
int* get_score_cell(int h[][max_adapter_size], int i, int j, std::string read,
                    std::string adap) {
  int s[] = {0, 0, 0};
  int* tmp = new int[2];
  tmp[0] = 0;
  tmp[1] = 0;
  s[0] = h[i - 1][j - 1] + get_sub_score(read[i - 1], adap[j - 1]);
  s[1] = h[i - 1][j] - get_gap_penalty(1, adap[j - 1]);
  s[2] = h[i][j - 1] - get_gap_penalty(1, read[i - 1]);
  int max = (s[0] >= s[1]) ? s[0] : s[1];
  max = (max >= s[2]) ? max : s[2];
  max = (max > 0) ? max : 0;
  tmp[0] = max;
  if (max == s[0]) {
    tmp[1] = 1;
  } else if (max == s[1]) {
    tmp[1] = 3;
  } else if (max == s[2]) {
    tmp[1] = 2;
  }
  return tmp;
}

void print_alignment(char a[2][max_read_size], int n) {
  std::cout << "Alignment: " << std::endl << std::endl;
  for (int j = 0; j < 2; j++) {
    if (j == 0)
      std::cout << "Read: ";
    else if (j == 1)
      std::cout << "Adap: ";
    for (int i = n - 1; i >= 0; i--) {
      std::cout << a[j][i] << " ";
    }
    std::cout << std::endl;
  }
}

int* align_seqs(std::string read, std::string adap) {
  int h[max_read_size][max_adapter_size], t[max_read_size][max_adapter_size],
      m = read.length() + 1, n = adap.length() + 1, max_i = 0, max_j = 0,
      max_v = -1, max_score, *rt = new int[2], *tmp = new int[2];
  rt[0] = -1;
  rt[1] = read.length();
  // Initialize first column to 0
  for (int i = 0; i < n; i++) h[0][i] = 0;
  // Initialize first row to 0
  for (int i = 0; i < m; i++) h[i][0] = 0;
  // print_matrix(h, m, n, read, adap);
  for (int i = 1; i < m; i++) {
    for (int j = 1; j < n; j++) {
      tmp = get_score_cell(h, i, j, read, adap);
      h[i][j] = *tmp;
      t[i][j] = *(tmp + 1);
      if (*tmp >= max_v) {
        max_i = i;
        max_j = j;
        max_v = *tmp;
      }
    }
  }
  max_score = max_v;
  // if(max_i != read.length())	// Ensure alignment to end of read.
  //   return rt;
  // std::cout << "Score: " << max_v << " ";
  // print_matrix(h, m, n, read, adap);
  // print_matrix(t, m, n, read, adap);
  // Traceback
  int _t, _l, _d, _align_n = 0;
  char _align[2][max_read_size];
  // 1 - Diag, 2->Right, 3->Down
  while (max_v != 0) {
    _d = h[max_i - 1][max_j - 1];
    _l = h[max_i][max_j - 1];
    _t = h[max_i - 1][max_j];
    switch (t[max_i][max_j]) {
      case 1:
        max_v = _d;
        max_i = max_i - 1;
        max_j = max_j - 1;
        _align[0][_align_n] = read[max_i];
        _align[1][_align_n] = adap[max_j];
        _align_n++;
        break;
      case 2:
        max_v = _l;
        max_j = max_j - 1;
        _align[0][_align_n] = '-';
        _align[1][_align_n] = adap[max_j];
        _align_n++;
        break;
      case 3:
        max_v = _t;
        max_i = max_i - 1;
        _align[0][_align_n] = read[max_i];
        _align[1][_align_n] = '-';
        _align_n++;
        break;
      case 0:
        max_v = 0;
        break;
    }
  }
  rt[1] = max_i;
  print_alignment(_align, _align_n);
  if (max_score <= ((_align_n - 2) * unit_score) - get_gap_penalty(2, 'A') ||
      max_j != 0)  // If alignment does not start at beginning of adapter
    return rt;
  rt[0] = max_score;
  // print_alignment(_align, _align_n);
  return rt;
}

// int find_adapters_contaminants(std::istream &cin, std::string
// adp_cntms_file){
//   std::vector<std::string> adp = read_adapters_from_fasta(adp_cntms_file,
//   "Read1"); std::string line, bases, qual, b; int ctr = 0, score, *t = new
//   int[2]; while (std::getline(cin, line)){
//     if(ctr % 400000 == 0)
//       std::cout << "Processed " << (ctr/4) << " reads ..." << std::endl;
//     if(ctr % 2 == 0){
//       ctr++;
//       continue;
//     }
//     if((ctr+1) % 2 == 0 && (ctr + 1) % 4 == 0 ){
//       qual = line;
//       // std::cout << bases << std::endl;
//       // std::cout << qual << std::endl;
//       score = -1;
//       for(std::vector<std::string>::iterator it = adp.begin(); it !=
//       adp.end() && score == -1; ++it) {
// 	b = ((*it).length() > bases.length()) ? bases :
// bases.substr(bases.length() - (*it).length(), (*it).length()); 	t =
// align_seqs(b, *it); 	score = *t; 	if(score != -1) 	  std::cout << "Trimmed: " <<
// bases.substr(0, bases.length() - *(t+1)) << std::endl;
//       }
//     } else {
//       bases = line;
//     }
//     ctr++;
//   }
//   return 0;
// }

// int main(int argc, char* argv[]){
//   std::string adp = "../data/adapters/NexteraPR.fa";
//   // clock_t begin = clock();
//   find_adapters_contaminants(std::cin, adp);
//   // clock_t end = clock();
//   // double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//   // std::cout << "Time Taken: " << elapsed_secs << std::endl;
//   return 0;
// }
