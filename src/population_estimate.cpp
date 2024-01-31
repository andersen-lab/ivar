#include "population_estimate.h"
#include <fstream>
#include <iostream>
#include <sstream>

float convert_decimal_date(std::string date){
  float decimal_date=0;
  //std::vector<float> parse_string;
  std::istringstream ss(date);

  std::string tmp;
  while(std::getline(ss, tmp, '-')){
    //words.push_back(tmp);
    std::cerr << tmp << std::endl;
  }
  return(decimal_date);
}

int estimate_populations(std::string variants, float evol_rate, std::string sample_date, std::string ref_date){
  int retval = 0;
  std::cerr << evol_rate << " " << variants << " " << sample_date << " " << ref_date << std::endl; 
  float sample_decimal = convert_decimal_date(sample_date); 
  std::cerr << sample_decimal << std::endl;
  /*
  uint32_t num_var = 0;
  for(uint32_t i=0; i < variants.size(); i++){
    if(!variants[i].amplicon_flux && !variants[i].depth_flag && variants[i].freq > 0.01 && !variants[i].is_ref){
      num_var += 1;
    }
  }
  float percent_possible_mutations = (num_var * beta) / (ref * 3);
  float expected_nt_sub_per_site = evolutionary_rate * time_elapsed;
  float num_populations = percent_possible_mutations / expected_nt_sub_per_site;

  //can't have 0 or fewer populations present
  uint32_t n = round(num_populations);
  if (n-1 <= 0) {
    population_bounds.push_back(1);
  } else {
    population_bounds.push_back(n-1);
  }
  population_bounds.push_back(n+1);
  */
  return(retval);
}
