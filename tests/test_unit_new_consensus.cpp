#include <iostream>
#include <vector>
#include <fstream>
#include "../src/include/armadillo"
#include "htslib/sam.h"
#include "../src/gmm.h"
#include "../src/saga.h"
#include "../src/ref_seq.h"
#include "../src/parse_gff.h"
#include "../src/call_consensus_clustering.h"
#include "../src/estimate_error.h"
#include "../src/solve_clustering.h"
#include "../src/interval_tree.h"

int main() {
  int num_tests = 4;
  int success = 0;

  /* TEST 1 - Position level masking based on amplicon fluctuation.
   * If a position is on an amplicon experiencing fluctuation, \
   * and all amplicon-specific variant frequencies are not assigned \
   * to the same cluster, we keep the position marked as experiencing \
   * fluctuation.
   *
   * First, we test to make sure amplicon-specific variant frequencies are \
   * assigned to the correct clusters.
   * Second, we check the amplicon flagging.
   *
   * A - Position not flagged.
   * B - Position flagged.
   */
  std::vector<variant> variants;
  gaussian_mixture_model model;
  arma::gmm_diag gmm_model;  //create the actual armadillo model

  std::vector<double> means = {0.90, 0.10};
  arma::mat mmeans(means);
  mmeans = mmeans.t();

  arma::mat covs = { {0.001, 0.001} };
  arma::mat hefts = { {0.5, 0.5} };
  gmm_model.set_params(mmeans, covs, hefts);
  model.model = gmm_model;
  model.n = 2;
  model.means = means;

  variant tmp;
  tmp.freq = 0.90;
  tmp.amplicon_flux = true;
  tmp.freq_numbers = {0.80, 1};
  variants.push_back(tmp);
  amplicon_specific_cluster_assignment(variants, model);

  variant output = variants[0];
  //this checks for correct assignment on the amplicon level
  if(output.freq_assignments[0] == 0 && output.freq_assignments[1] == 0){
    success++;
  }
  rewrite_position_masking(variants); //change the position masking
  if(!variants[0].amplicon_flux) success++; //should not be flagged


  variants.clear();
  means.clear();
  means = {0.60, 0.40};
  mmeans = means;
  mmeans = mmeans.t();
  gmm_model.set_params(mmeans, covs, hefts);
  model.model = gmm_model;
  model.means = means;
  tmp.freq = 0.5;
  tmp.amplicon_flux = true;
  tmp.freq_numbers = {0.55, 0.45};
  variants.push_back(tmp);
  amplicon_specific_cluster_assignment(variants, model);

  output = variants[0];
  //here they should be assigned to different clusters
  if(output.freq_assignments[0] == 0 && output.freq_assignments[1] == 1){
    success++;
  }
  if(variants[0].amplicon_flux){
    success++;
  }

  /* TEST 2 - Amplicons are masked only if they
   * (a) they contain a position flagged as experiencing fluctuation
   * (b) a cluster exists that is within the standard deviation bounds of the cluster the variant is assigned to
   * In other words if we know a position experiences a std devation \
   * of 0.10 and it's assigned to a cluster at 0.50, however a cluster \
   * also exists at 0.45 we would flag the amplicon.
   */

  //std::vector<uint32_t> amplicons_to_mask = rewrite_amplicon_masking(variants, means);



  std::cerr << "success " << success << " tests " << num_tests << std::endl;
  return (num_tests == success) ? 0 : -1;
}
