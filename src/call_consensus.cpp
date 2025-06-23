#include <stdint.h>

#include <fstream>
#include <iostream>

#include "htslib/vcf.h"

struct maj {
  int *ind;
  int n;
};

#define GET_MAX_INDICES(type_t, pu, no, m) \
  {                                        \
    type_t v = 1;                          \
    type_t *p = (type_t *)(pu);            \
    for (int k = 0; k < no; ++k) {         \
      if (p[k] > v) {                      \
        v = p[k];                          \
        m.n = 0;                           \
        m.ind[m.n] = k;                    \
      } else if (p[k] == v) {              \
        m.n++;                             \
        m.ind[m.n] = k;                    \
      }                                    \
    }                                      \
  }

#define INDEL_SIZE 30

#define GAP 'N'

// From bcftools.h -
// https://github.com/samtools/bcftools/blob/b0376dff1ed70603c9490802f37883b9009215d2/bcftools.h#L48
static inline char gt2iupac(char a, char b) {
  static const char iupac[4][4] = {{'A', 'M', 'R', 'W'},
                                   {'M', 'C', 'S', 'Y'},
                                   {'R', 'S', 'G', 'K'},
                                   {'W', 'Y', 'K', 'T'}};
  if (a >= 'a') a -= 'a' - 'A';
  if (b >= 'a') b -= 'a' - 'A';
  if (a == 'A')
    a = 0;
  else if (a == 'C')
    a = 1;
  else if (a == 'G')
    a = 2;
  else if (a == 'T')
    a = 3;
  else
    return 'N';
  if (b == 'A')
    b = 0;
  else if (b == 'C')
    b = 1;
  else if (b == 'G')
    b = 2;
  else if (b == 'T')
    b = 3;
  else
    return 'N';
  return iupac[(int)a][(int)b];
}

void get_consensus_indel(char **allele, maj m, char *&nuc) {
  int max_len = 0, l = 0;
  uint8_t q = 0;
  for (int k = 0; k <= m.n; ++k) {
    if (strlen(allele[m.ind[k]]) > max_len) max_len = strlen(allele[m.ind[k]]);
  }
  char n;
  for (int j = 0; j < max_len; ++j) {  // Iterate over base positions
    n = 0;
    for (int k = 0; k <= m.n; ++k) {  // Iterate over major alleles
      if (j < strlen(allele[m.ind[k]])) {
        if (n == 0) {
          n = allele[m.ind[k]][j];
        } else {
          n = gt2iupac(n, allele[m.ind[k]][j]);
        }
      }
    }
    if (n == 0)
      nuc[l] = GAP;
    else
      nuc[l] = n;
    l++;
    if (l % INDEL_SIZE == 0) {
      nuc = (char *)realloc(nuc,
                            INDEL_SIZE * ((l / INDEL_SIZE) + 1) * sizeof(char));
    }
  }
  nuc[l] = 0;
}

int main_old(int argc, char *argv[]) {
  std::cout << "Path " << argv[1] << std::endl;
  std::string bcf = std::string(argv[1]);
  std::string region_;
  std::string out_path = std::string(argv[2]);
  if (argc > 3) {
    region_ = std::string(argv[3]);
  }
  if (!bcf.empty()) {
    // open BCF for reading
    htsFile *in = bcf_open(bcf.c_str(), "r");
    std::ofstream out(out_path.c_str());
    out << ">" << out_path << "\n";
    if (in == NULL) {
      throw std::runtime_error("Unable to open BCF/VCF file.");
    }
    int ad_id;
    // Get the header
    bcf_hdr_t *header = bcf_hdr_read(in);
    if (header == NULL) {
      bcf_close(in);
      throw std::runtime_error("Unable to open BCF header.");
    }
    ad_id = bcf_hdr_id2int(header, BCF_DT_ID, "AD");
    std::cout << ad_id << " " << std::endl;
    // Initialize iterator
    bcf1_t *v = bcf_init();
    int ctr = 0, prev_pos = 0, s, d;
    maj m;
    char *nuc = (char *)malloc(INDEL_SIZE * sizeof(char)), *r;
    while (bcf_read(in, header, v) == 0) {
      // std::cout << "Position: " << v->pos << std::endl;
      // std::cout << "Number of Alleles: " << v->n_allele << std::endl;
      for (int i = prev_pos; i < v->pos - 1; ++i) {
        out << GAP;
      }
      bcf_unpack(v, BCF_UN_FMT);
      bcf_unpack(v, BCF_UN_STR);
      bcf_fmt_t *fmt = v->d.fmt;
      for (int i = 0; i < v->n_fmt; ++i) {
        if (strcmp(header->id[BCF_DT_ID][fmt[i].id].key, "AD") == 0) {
          m.ind = (int *)malloc(fmt[i].n * sizeof(int));
          m.n = -1;
          switch (fmt[i].type) {
            case BCF_BT_INT8:
              GET_MAX_INDICES(int8_t, fmt[i].p, fmt[i].n, m);
              break;
            case BCF_BT_INT16:
              GET_MAX_INDICES(int16_t, fmt[i].p, fmt[i].n, m);
              break;
            case BCF_BT_INT32:
              GET_MAX_INDICES(int32_t, fmt[i].p, fmt[i].n, m);
              break;
            case BCF_BT_FLOAT:
              GET_MAX_INDICES(float, fmt[i].p, fmt[i].n, m);
              break;
            default:
              hts_log_error("Unexpected type %d", fmt[i].type);
              exit(1);
          }
          if (m.n == -1) {
            out << GAP;
          } else {
            if (bcf_get_variant_types(v) == VCF_REF) {
              out << v->d.allele[0];
            } else if (bcf_get_variant_types(v) == VCF_SNP ||
                       bcf_get_variant_types(v) == VCF_MNP) {
              nuc = (char *)malloc(INDEL_SIZE * sizeof(char));
              get_consensus_indel(v->d.allele, m, nuc);
              out << nuc;
            } else {  // For INDEL insert difference.between reference and
                      // consensus allele
              nuc = (char *)malloc(INDEL_SIZE * sizeof(char));
              get_consensus_indel(v->d.allele, m, nuc);
              s = strlen(v->d.allele[0]) -
                  strlen(nuc);  //  > 0 Deletion. < 0 INSERTION
              d = (s > 0) ? strlen(v->d.allele[0]) : strlen(nuc);
              r = (s > 0) ? v->d.allele[0] : nuc;
              s = (s < 0) ? -1 * s : s;
              if (v->pos == 499)
                std::cout << "Number: " << m.n
                          << " Maj Allele: " << v->d.allele[m.ind[0]]
                          << strlen(nuc) << " - Nuc Length" << nuc << " " << s
                          << " " << d << std::endl;
              for (int i = 0; i < s; ++i) {
                out << r[d - s];
              }
              free(m.ind);
              free(nuc);
              // std::cout << nuc << std::endl;
            }
          }
        }
      }
      prev_pos = v->pos;
      ctr++;
      if (ctr % 1000 == 0) {
        std::cout << ctr << std::endl;
      }
    }
    bcf_destroy(v);
    bcf_hdr_destroy(header);
    bcf_close(in);
    out.close();
  }
  return 0;
}
