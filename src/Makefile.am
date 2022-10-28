LIBS = -lhts -lz -lpthread

CXXFLAGS = -v -g -std=c++11 -I./alglib -Wall -Wextra -Werror -Wno-unused-parameter -Wno-unneeded-internal-declaration

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = ivar
ivar_SOURCES = ivar.cpp call_consensus_pileup.cpp alignment.cpp suffix_tree.cpp trim_primer_quality.cpp remove_reads_from_amplicon.cpp call_variants.cpp primer_bed.cpp allele_functions.cpp get_masked_amplicons.cpp get_common_variants.cpp parse_gff.cpp ref_seq.cpp interval_tree.cpp clustering.cpp kmeans.cpp ./alglib/ap.cpp ./alglib/dataanalysis.cpp ./alglib/alglibinternal.cpp ./alglib/alglibmisc.cpp ./alglib/diffequations.cpp ./alglib/fasttransforms.cpp ./alglib/integration.cpp ./alglib/interpolation.cpp ./alglib/kernels_avx2.cpp ./alglib/kernels_fma.cpp ./alglib/kernels_sse2.cpp ./alglib/linalg.cpp ./alglib/optimization.cpp ./alglib/solvers.cpp ./alglib/specialfunctions.cpp ./alglib/statistics.cpp 
