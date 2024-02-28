#include "population_estimate.h"
#include "gmm.h"
#include <numeric>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>

float convert_decimal_date(std::string date){
  std::istringstream date_stream(date);
  std::tm tm = {};
  date_stream >> std::get_time(&tm, "%m-%d-%Y");
    
  if (date_stream.fail()) {
    std::cerr << "Error parsing the date format." << std::endl;
    return 0.0;
  }

  //convert the date to a decimal value
  int year = tm.tm_year + 1900;
  float decimal_date = year + (tm.tm_yday + 1) / 365.25;
  return(decimal_date);
}

void linear_regression(const std::vector<float>& x_values, const std::vector<float>& y_values,
                       float& slope, float& intercept) {
    size_t n = x_values.size();

    // Calculate the means of x and y
    float mean_x = std::accumulate(x_values.begin(), x_values.end(), 0.0) / n;
    float mean_y = std::accumulate(y_values.begin(), y_values.end(), 0.0) / n;

    // Calculate the slope (m) and intercept (b)
    float numerator = 0.0, denominator = 0.0;
    for (size_t i = 0; i < n; ++i) {
        numerator += (x_values[i] - mean_x) * (y_values[i] - mean_y);
        denominator += (x_values[i] - mean_x) * (x_values[i] - mean_x);
    }
    slope = numerator / denominator;
    intercept = mean_y - slope * mean_x;
}

int estimate_populations(std::string variants, float evol_rate, std::string sample_date, std::string ref_date){
  int retval = 0;
  float lower_bound = 0.03;
  float upper_bound = 0.97;
  uint32_t depth_cutoff = 10;
  float quality_threshold = 20;
  uint32_t round_val = 3;
  std::vector<variant> variant_vec;
  float sample_decimal = convert_decimal_date(sample_date); 
  float reference_decimal = convert_decimal_date(ref_date);
  if(reference_decimal > sample_decimal){
    std::cerr << "Error reference date occurs after sample date." << std::endl;
    return(-1);
  }
  float time_elapsed = sample_decimal - reference_decimal;
  std::cerr << time_elapsed << std::endl;
  std::vector<uint32_t> low_quality_positions = find_low_quality_positions(variants, depth_cutoff, lower_bound, upper_bound, quality_threshold, round_val);
  std::vector<uint32_t> deletion_positions = find_deletion_positions(variants, depth_cutoff, lower_bound, upper_bound, round_val);
  parse_internal_variants(variants, variant_vec, depth_cutoff, lower_bound, upper_bound, deletion_positions, low_quality_positions, round_val); 

  std::vector<variant> kept_variants;
  float var_count = 0;
  float ref_length = 0;
  for(uint32_t i=0; i < variant_vec.size(); i++){
    if((float)variant_vec[i].position > ref_length){
      ref_length = (float)variant_vec[i].position;
    }
    if(!variant_vec[i].amplicon_flux && !variant_vec[i].depth_flag && !variant_vec[i].is_ref && variant_vec[i].freq > lower_bound && variant_vec[i].qual > 20){
      std::cerr << variant_vec[i].nuc << " " << variant_vec[i].freq << " " << variant_vec[i].position << " " << variant_vec[i].qual << std::endl;
      kept_variants.push_back(variant_vec[i]);
      var_count += 1;
    }
  }
  std::vector<float> x;
  std::vector<float> y; 
  //regress
  for(uint32_t i = 2; i < 6; i++){
    float section = time_elapsed / float(i);
    float cumulative_mutations = 0;
    for(uint32_t j=1; j < i+1; j++){
      float time_point = section * (float)j;
      float mut_expected = (time_point * evol_rate * ref_length);
      cumulative_mutations += mut_expected;
    }
    x.push_back((float)i);
    y.push_back(cumulative_mutations);
  }
  float slope, intercept;

  // Calculate linear regression
  linear_regression(x, y, slope, intercept);
  float pop = (var_count-intercept) / slope;
  std::cerr << pop << std::endl;
  return(retval);
}
