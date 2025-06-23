#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#ifndef primer_bed
#define primer_bed

/* Substitution Matrix

   A  T  G  C N
A  1 -1 -1 -1 0
T -1  1 -1 -1 0
G -1 -1  1 -1 0
C -1 -1 -1  1 0
N  0  0  0  0 0
*/

const int unit_score = 2;

const int substitution[5][5] = {
    {unit_score, -unit_score, -unit_score, -unit_score, 0},
    {-unit_score, unit_score, -unit_score, -unit_score, 0},
    {-unit_score, -unit_score, unit_score, -unit_score, 0},
    {-unit_score, -unit_score, -unit_score, unit_score, 0},
    {0, 0, 0, 0, 0}};

const int gap_open = unit_score - 1;
const int gap_extension = -1;
const int max_read_size = 500;
const int max_adapter_size = 60;
const int MAX_GAP = 2;

int get_sub_score(char a, char b);
int get_gap_penalty(int k, char a);
void print_matrix(int h[][max_adapter_size], int r, int c, std::string read,
                  std::string adap);
int* get_score_cell(int h[][max_adapter_size], int i, int j, std::string read,
                    std::string adap);
void print_alignment(char a[2][max_read_size], int n);
int* align_seqs(std::string read, std::string adap);
int find_adapters_contaminants(std::istream& cin, std::string adp_cntms_file);
int main(int argc, char* argv[]);

#endif
