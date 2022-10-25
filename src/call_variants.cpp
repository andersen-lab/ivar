#include "call_variants.h"

const float sig_level = 0.01;

std::vector<allele>::iterator get_ref_allele(std::vector<allele> &ad, char ref){
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc[0] == ref)
      return it;
  }
  // print_allele_depths(ad);
  return ad.end();
}

double* get_frequency_depth(allele a, uint32_t pos_depth, uint32_t total_depth){ // pos_depth is ungapped depth after passing quality. total_depth is depth without quality filtering for indels.
  double *val = new double[2];
  if(a.nuc[0] == '+' || a.nuc[0] == '-'){	// For insertions and deletions use depth discarding quality
    val[0] = a.depth/(double)total_depth;
    val[1] = total_depth;
    return val;
  }
  val[0] = a.depth/(double)pos_depth;
  val[1] = pos_depth;
  return val;
}

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double min_threshold, uint8_t min_depth, std::string ref_path, std::string gff_path){
  std::string line, cell, bases, qualities, region;
  ref_antd refantd(ref_path, gff_path);
  std::ostringstream out_str;
  std::ofstream fout((out_file+".tsv").c_str());
  //make columns in the file
  fout << "REGION"
    "\tPOS"
    "\tREF"
    "\tALT"
    "\tREF_DP"
    "\tREF_RV"
    "\tREF_QUAL"
    "\tALT_DP"
    "\tALT_RV"
    "\tALT_QUAL"
    "\tALT_FREQ"
    "\tTOTAL_DP"
    "\tPVAL"
    "\tPASS"
    "\tGFF_FEATURE"
    "\tREF_CODON"
    "\tREF_AA"
    "\tALT_CODON"
    "\tALT_AA"
    "\tPOS_AA"
       << std::endl;
  int ctr = 0;
  int64_t pos = 0;
  uint32_t mdepth = 0, pdepth = 0, gdepth =0; // mpdepth for mpileup depth and pdeth for ungapped depth at position
  double pval_left, pval_right, pval_twotailed, *freq_depth, err;
  std::stringstream line_stream;
  char ref;
  std::vector<allele> ad;
  std::vector<allele>::iterator ref_it;
  while (std::getline(cin, line)){
    line_stream << line;
    ctr = 0;
    ref = 'N';
    while(std::getline(line_stream,cell,'\t')){
      switch(ctr){
      case 0:
	region = cell;
	break;
      case 1:
	pos = stoi(cell);
	break;
      case 2:
	ref = refantd.get_base(pos, region);
	ref = (ref == 0) ? cell[0] : ref; // If ref does not exist then use from mpileup
	break;
      case 3:
	mdepth = stoi(cell);
	break;
      case 4:
	bases = cell;
	break;
      case 5:
	qualities = cell;
	break;
      case 6:
	break;
      }
      ctr++;
    }
    
    ad = update_allele_depth(ref, bases, qualities, min_qual);
    if(ad.size() == 0){
      line_stream.clear();
      continue;
    }
    // Get ungapped depth
    pdepth = 0;
    //get gapped depth
    gdepth = 0;
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      gdepth += it->depth;
      //commenting these two lines gets the gapped depth
      if(it->nuc[0]=='*' || it->nuc[0] == '+' || it->nuc[0] == '-')
	continue;
      pdepth += it->depth;
    }
    if(pdepth < min_depth) {	// Check for minimum depth
      line_stream.clear();
      continue;
    }
    ref_it = get_ref_allele(ad, ref);
    if (ref_it == ad.end()) {	// If ref not present in reads.
      allele a;
      a.nuc = ref;
      a.depth = 0;
      a.reverse = 0;
      a.mean_qual = 0;
      ad.push_back(a);
      ref_it = ad.end() - 1;
    }
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if((*it == *ref_it) || it->nuc[0]=='*')
	continue;
      //use gapped depth if we're looking at insertion or deletion
      if (it->nuc[0] == '+' || it->nuc[0] == '-'){
        freq_depth = get_frequency_depth(*it, gdepth, gdepth);
      }else{
        freq_depth = get_frequency_depth(*it, pdepth, mdepth);
      }
      if(freq_depth[0] < min_threshold)

	continue;
      //this only adds the first bit of the information to the tsv
      out_str << region << "\t";
      out_str << pos << "\t";
      out_str << ref << "\t";
      out_str << it->nuc << "\t";
      out_str << ref_it->depth << "\t";
      out_str << ref_it->reverse << "\t";
      out_str << (uint16_t)ref_it->mean_qual << "\t";
      out_str << it->depth << "\t";
      out_str << it->reverse << "\t";
      out_str << (uint16_t) it->mean_qual << "\t";
      out_str << freq_depth[0] << "\t";
      out_str << freq_depth[1] << "\t";
      /*
	    | Var   | Ref      |
	Exp | Error | Err free |
	Obs | AD    | RD       |
       */
      //adding pvalue related info
      err = pow(10, ( -1 * (it->mean_qual)/10));
      kt_fisher_exact((err * mdepth), (1-err) * mdepth, it->depth, ref_it->depth, &pval_left, &pval_right, &pval_twotailed);
      out_str << pval_left << "\t";
      if(pval_left <= sig_level){
	out_str << "TRUE" << "\t";
      } else {
	out_str << "FALSE" << "\t";
      }
      //adding codon related info, can add gene info here too
      if(it->nuc[0] != '+' && it->nuc[0] != '-'){
	refantd.codon_aa_stream(region, out_str, fout, pos, it->nuc[0]);
      } else {
	fout << out_str.str() << "NA\tNA\tNA\tNA\tNA\tNA" << std::endl;
      }
      out_str.str("");
      out_str.clear();
      delete[] freq_depth;
    }
    line_stream.clear();
  }
  fout.close();
  return 0;
}
