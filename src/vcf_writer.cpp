#include "vcf_writer.h"

int vcf_writer::init_header() {
  this->hdr = bcf_hdr_init("w");
  int res;
  res = bcf_hdr_append(
      hdr,
      "##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
  res = bcf_hdr_append(
      hdr,
      "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
  res = bcf_hdr_append(
      hdr,
      "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for "
      "each allele (based on minimum quality threshold)\">");
  res = bcf_hdr_append(
      hdr,
      "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total read depth for "
      "each allele on forward strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(
      hdr,
      "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total read depth for "
      "each allele on reverse strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(
      hdr,
      "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for "
      "each ALT allele in the same order as listed (estimated from primary "
      "data, not called genotypes)\">");
  res = bcf_hdr_append(
      hdr, ("##contig=<ID=cref,assembly=" + this->region + ">").c_str());
  res = bcf_hdr_append(
      hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  res = bcf_hdr_append(hdr,
                       "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="
                       "\"Genotype Quality\">");
  res = bcf_hdr_append(hdr,
                       "##FILTER=<ID=fobs,Description=\"Fisher's exact test to "
                       "check if frequency of the iSNV is significantly higher "
                       "than the mean error rate at that position\">");
  bcf_hdr_add_sample(hdr, this->sample_name.c_str());
  if (res != 0) return -1;
  return 0;
}

int vcf_writer::init(char _mode, std::string fname, std::string region,
                     std::string sample_name, std::string ref_path) {
  this->ref = new ref_antd(ref_path);
  this->region = region;
  this->sample_name = sample_name;
  std::string mode = "w";
  mode += _mode;
  vcfFile *vcf_file = bcf_open(fname.c_str(), mode.c_str());
  this->init_header();
  if (hdr == NULL || vcf_file == NULL) {
    std::cout << "Unable to write BCF/VCF file" << std::endl;
    return -1;
  }
  if (bcf_hdr_write(vcf_file, hdr) < 0) {
    std::cout << "Unable to write BCF/VCF file" << std::endl;
    return -1;
  }
  return 0;
}

int vcf_writer::write_record(uint32_t pos, std::vector<allele> &alleles,
                             char ref) {
  bcf1_t *rec = bcf_init();
  rec->rid = rec->rid = bcf_hdr_name2id(this->hdr, this->region.c_str());
  rec->pos = pos;
  allele aref = alleles.at(find_ref_in_allele(alleles, ref));
  rec->qual = aref.mean_qual;
}

int main(int argc, char *argv[]) {
  vcfFile *vcf_file = bcf_open("../data/test.vcf", "w");
  bcf_hdr_t *h = bcf_hdr_init("w");
  int res;

  if (res < 0) std::cout << "Unable to write BCF header" << std::endl;
  // Add sample
  bcf_hdr_add_sample(h, "TEST");
  // bcf_hdr_add_sample(h, NULL);
  bcf_hdr_write(vcf_file, h);
  // Write record
  bcf1_t *rec = bcf_init();
  rec->rid = rec->rid = bcf_hdr_name2id(h, "20");
  rec->pos = 0;
  rec->qual = 20;
  // ID
  bcf_update_id(h, rec, "rs6054257");
  // INFO
  bcf_update_alleles_str(h, rec, "G,A,TA");
  // FILTER
  int32_t tmpi = bcf_hdr_id2int(h, BCF_DT_ID, "PASS");
  bcf_update_filter(h, rec, &tmpi, 1);
  // FORMAT
  int32_t *tmpia = (int *)malloc(bcf_hdr_nsamples(h) * 2 * sizeof(int));
  tmpia[0] = bcf_gt_phased(0);
  tmpia[1] = bcf_gt_phased(0);
  bcf_update_genotypes(h, rec, &tmpi, bcf_hdr_nsamples(h));
  tmpia[0] = 12;
  tmpia[1] = 4;
  bcf_update_format_int32(h, rec, "GQ", tmpia, bcf_hdr_nsamples(h));
  tmpi = 1;
  bcf_update_info_int32(h, rec, "NS", &tmpi, 1);
  if (res < 0) std::cout << "Alleles not updated" << std::endl;
  int32_t td = 12;
  int32_t ad = 4;
  bcf_update_info_int32(h, rec, "DP", &td, 1);
  bcf_update_info_int32(h, rec, "AD", &ad, 1);
  std::cout << bcf_hdr_nsamples(h) << std::endl;
  std::cout << rec->n_sample << std::endl;
  bcf_write1(vcf_file, h, rec);
  bcf_hdr_destroy(h);
  bcf_close(vcf_file);
  return 0;
}
