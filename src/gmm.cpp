#include "./include/armadillo"
#include <fstream>
#include <cmath>

void go(uint32_t offset, uint32_t k, std::vector<double> means, std::vector<double> combination, std::vector<std::vector<double>> &combos, double error) {
  if (k == 0) {
    double sum_of_elems = std::accumulate(combination.begin(), combination.end(), 0.0);
    std::cerr << sum_of_elems << std::endl;
    if(sum_of_elems < 1+error && sum_of_elems > 1-error){
      combos.push_back(combination);
    }
    return;
  }
  for (uint32_t i = offset; i <= means.size() - k; ++i) {
    combination.push_back(means[i]);
    go(i+1, k-1, means, combination, combos, error);
    combination.pop_back();
  }
}
void solve_solution_sets(std::vector<double> means, uint32_t n){
  if(means.size() < n){
    return;
  }
  //generate the combinations
  std::vector<double> combination; //placeholder to use in function
  std::vector<std::vector<double>> combos;
  double errors = 0.05
  go(0, n, means, combination, combos, error);
  for(uint32_t i=0; i < combos.size(); i++){
    for(uint32_t j=0; j < combos[i].size(); j++){
      std::cerr << combos[i][j] << " ";
    }
    std::cerr << std::endl;
  }

}

void determine_outlier_variants(std::vector<std::vector<double>> prob_matrix, std::vector<float> filtered_frequencies, std::vector<uint32_t> &outlier_idxs){
  /*
   * Here we take in prob of every variants being assigned to every gaussian and determine which variants can be concretely assigned and which cannot.
   */
  double threshold = 3;
  //iterate the transposed version
  for(uint32_t i=0; i < prob_matrix[0].size(); i++) {
    std::vector<double> tmp;
    for(uint32_t j=0; j < prob_matrix.size(); j++){
      tmp.push_back(exp(prob_matrix[j][i]));
    }
    std::sort(tmp.begin(),tmp.end(), std::greater<double>());
    double ratio = tmp[0] / tmp[1];
    if (tmp[1] * threshold > tmp[0]){
      outlier_idxs.push_back(i);
      std::cerr << ratio << " " << filtered_frequencies[i] << std::endl;
    }
  }
}

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
  //the n min and n max represent the range of n values
  uint32_t n = 2;

  //filter out the data we won't use in our initial model
  std::vector<float> filtered_frequencies;
  for(uint32_t i=0; i < frequencies.size(); i++){
    if(depths[i] > 10 && flagged[i] == "FALSE" && frequencies[i] > lower_bound && frequencies[i] < upper_bound){
      filtered_frequencies.push_back(frequencies[i]);  
    }
  }

  //initialize armadillo dataset and populate with frequency data
  //(rows, cols) where each columns is a sample
  arma::mat data(1, filtered_frequencies.size(), arma::fill::zeros); 
  for(uint32_t i = 0; i < filtered_frequencies.size(); i++){
    double tmp = static_cast<double>(filtered_frequencies[i]);
    data.col(i) = tmp;
  }
  //model learning
  bool status = model.learn(data, n, arma::eucl_dist, arma::random_subset, 10, 5, 1e-10, true);
  if(status == false){
    std::cerr << "gmm model failed" << std::endl;
  }
  //get the means of the gaussians
  std::vector<double> means;
  
  //get the probability of each frequency being assigned to each gaussian
  std::vector<std::vector<double>> prob_matrix;
  for(uint32_t i=0; i < n; i++){
    means.push_back((double)model.means[i]);
    arma::rowvec set_likelihood = model.log_p(data.cols(0,filtered_frequencies.size()-1), i);
    std::vector<double> tmp;
    for(uint32_t j=0; j < filtered_frequencies.size(); j++){
      tmp.push_back((double)set_likelihood[j]);
    }
    prob_matrix.push_back(tmp);
  }
  
  //based on the probabilities, determine which points cannot be concretely assigned
  std::vector<uint32_t> outlier_idxs; 
  double raw_distance_threshold = 0.20; //here we define any point existing too far outside of a cluster as not possible
  determine_outlier_variants(prob_matrix, filtered_frequencies, outlier_idxs, raw_distance_threshold);
  
  //solve for solutions adding to ~1
  solve_solution_sets(means, n);
  
  //places where we'll call an N value
  //std::vector<uint32_t> unassigned_points;
  //model.save("my_model.gmm");
  return(retval);
}
