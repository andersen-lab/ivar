#include "./include/armadillo"
#include <fstream>


void split(const std::string &s, char delim, std::vector<std::string> &elems){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

void parse_internal_variants(std::string filename, std::vector<float> &frequencies, std::vector<uint32_t> &positions, std::vector<uint32_t> &depths, std::vector<std::string> &flagged, std::vector<std::string> &nucs){
  /*
   * Parses the variants file produced internally by reading bam file line-by-line.
   */
  std::ifstream infile(filename);
  std::string line;
  uint32_t count = 0;
  while (std::getline(infile, line)) {
    if(count == 0) {
      count += 1;
      continue;
    }
    std::vector<std::string> row_values;
    //split the line by delimiter
    uint32_t pos = 0;
    uint32_t depth = 0;
    float freq = 0;
    std::string flag = "";
    split(line, '\t', row_values);
    pos = std::stoi(row_values[0]);
    nucs.push_back(row_values[1]);
    depth = std::stoi(row_values[2]);
    freq = std::stof(row_values[3]);
    flag = row_values[5];
    frequencies.push_back(freq);
    positions.push_back(pos);
    depths.push_back(depth);
    flagged.push_back(flag);
    count += 1;
  } 
}

int gmm_model(std::string prefix){
  arma::gmm_diag model;
  int retval = 0;
  float lower_bound = 0.03;
  float upper_bound = 0.97;

  std::string filename = prefix + ".txt";
  std::vector<float> frequencies;
  std::vector<uint32_t> positions;
  std::vector<uint32_t> depths;
  std::vector<std::string> nucs;
  std::vector<std::string> flagged;
  parse_internal_variants(filename, frequencies, positions, depths, flagged, nucs);

  //filter out the data we won't use in our initial model
  std::vector<uint32_t> filtered_frequencies;
  for(uint32_t i=0; i < frequencies.size(); i++){
    if(depths[i] > 10 && flagged[i] == "FALSE" && frequencies[i] > lower_bound && frequencies[i] < upper_bound){
      filtered_frequencies.push_back(frequencies[i]);  
    }
  }
  //initialize armadillo dataset
  arma::fmat data(1, filtered_frequencies.size()); 

  return(retval);
}
