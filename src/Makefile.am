LIBS = -lhts -lz -lpthread

CXXFLAGS = -v -g -std=c++11 -Wall -Wextra -Werror

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = ivar
ivar_SOURCES = ivar.cpp call_consensus_pileup.cpp alignment.cpp suffix_tree.cpp trim_primer_quality.cpp remove_reads_from_amplicon.cpp call_variants.cpp primer_bed.cpp allele_functions.cpp get_masked_amplicons.cpp get_common_variants.cpp parse_gff.cpp ref_seq.cpp interval_tree.cpp saga.cpp include/armadillo_bits/*.hpp include/armadillo gmm.cpp population_estimate.cpp
 
ivar_LDADD = $(LIBS)
include_HEADERS = include/armadillo_bits/*.hpp
