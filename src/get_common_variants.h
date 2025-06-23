#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#ifndef get_common_variants
#define get_common_variants

int read_variant_file(
    std::ifstream &fin, unsigned int file_number,
    std::map<std::string, unsigned int> &counts,
    std::map<std::string, std::string> &file_tab_delimited_str);
int common_variants(std::string out, double min_threshold, char *files[],
                    unsigned int nfiles);

#endif
