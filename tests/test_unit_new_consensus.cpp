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
  int num_tests = 1;
  int success = 0;

  /* TEST 1 - Position level masking based on amplicon fluctuation.
   */
  std::vector<variant> variants;
  gaussian_mixture_model model;
  arma::gmm_diag gmm_model  //create the actual armadillo model
  model.model = gmm_model;
  model.n = 2;
  model.means = {0.90, 0.10};

  variant tmp;
  tmp.freq = 0.90;
  tmp.amplicon_flux = true;
  tmp.freq_numbers = {0.80, 1};
  variants.push_back(tmp);
  amplicon_specific_cluster_assignment(variants, model);

  variant output = variants[0];
  std::cerr << output.freq_assignments.size() << std::endl;
  exit(0);
  //rewrite_position_masking(std::vector<variant> &variants, gaussian_mixture_model model);

  /* TEST 2 - Amplicons are masked only if they
   * (a) they contain a position flagged as experiencing fluctuation
   * (b) a cluster exists that is within the standard deviation bounds of the cluster the variant is assigned to
   * In other words if we know a position experiences a std devation \
   * of 0.10 and it's assigned to a cluster at 0.50, however a cluster \
   * also exists at 0.45 we would flag the amplicon.
   */
  std::vector<double> means;
  //populate some fake variants to test
  std::vector<uint32_t> amplicons_to_mask = rewrite_amplicon_masking(variants, means);
  return (num_tests == success) ? 0 : -1;
}
