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
  int num_tests = 2;
  int success = 0;

  /* TEST 1 - Position level masking based on amplicon fluctuation.
   * If a position is on an amplicon experiencing fluctuation, \
   * and all amplicon-specific variant frequencies are not assigned \
   * to the same cluster, we keep the position marked as experiencing \
   * fluctuation.
   *
   * First, we test to make sure amplicon-specific variant frequencies are \
   * assigned to the correct clusters.
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
  if(output.freq_assignments[0] == 0 && output.freq_assignments[1] == 0){
    success++;
  }
  std::cerr << variants.size() << std::endl;
  rewrite_position_masking(variants); //change the position masking
  if(!variants[0].amplicon_flux) success++; //should not be flagged

  /* TEST 2 - Amplicons are masked only if they
   * (a) they contain a position flagged as experiencing fluctuation
   * (b) a cluster exists that is within the standard deviation bounds of the cluster the variant is assigned to
   * In other words if we know a position experiences a std devation \
   * of 0.10 and it's assigned to a cluster at 0.50, however a cluster \
   * also exists at 0.45 we would flag the amplicon.
   */
  //std::vector<uint32_t> amplicons_to_mask = rewrite_amplicon_masking(variants, means);
  return (num_tests == success) ? 0 : -1;
}
