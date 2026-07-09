add_executable(check_primer_trim test_primer_trim.cpp)

add_executable(check_trim test_trim.cpp)

add_executable(check_quality_trim test_quality_trim.cpp)

add_executable(check_consensus test_call_consensus_from_plup.cpp)

add_executable(check_allele_depth test_allele_depth.cpp)

add_executable(check_consensus_threshold test_consensus_threshold.cpp)

add_executable(check_consensus_min_insert_threshold test_consensus_min_insert_threshold.cpp)

add_executable(check_consensus_min_depth test_consensus_min_depth.cpp)

add_executable(check_consensus_seq_id test_consensus_seq_id.cpp)

add_executable(check_variants test_variants.cpp)

#add_executable(check_variants_2 test_new_variants.cpp)

add_executable(check_common_variants test_common_variants.cpp)

add_executable(check_primer_bed test_primer_bed.cpp)

add_executable(check_clustering_consensus test_clustering_consensus.cpp)

add_executable(check_clustering_majority test_clustering_majority.cpp)

add_executable(check_unit_consensus test_unit_new_consensus.cpp)

add_executable(check_getmasked test_getmasked.cpp)

add_executable(check_removereads test_removereads.cpp)

add_executable(check_unpaired_trim test_unpaired_trim.cpp)

add_executable(check_primer_trim_edge_cases test_primer_trim_edge_cases.cpp)

add_executable(check_isize_trim test_isize_trim.cpp)

add_executable(check_interval_tree test_interval_tree.cpp)

add_executable(check_amplicon_search test_amplicon_search.cpp)

add_executable(check_variant_caller test_variant_caller.cpp)

add_executable(check_site_aggregator test_site_aggregator.cpp)

add_executable(check_site_amplicon_aggregator_stats test_site_amplicon_aggregator_stats.cpp)

add_executable(check_site_state_aggregator_stats test_site_state_aggregator_stats.cpp)

add_executable(check_gmm_1d test_gmm_1d.cpp)

add_executable(check_gaussian_sampler test_gaussian_sampler.cpp)

add_executable(check_multinomial_sampler test_multinomial_sampler.cpp)

add_executable(check_bootstrap_variants test_bootstrap_variants.cpp)

add_executable(check_bd_threshold test_bd_threshold.cpp)

add_executable(check_subset_sum test_subset_sum.cpp)

add_executable(check_parse_gff test_parse_gff.cpp)

add_executable(check_ref_discontinuous_cds test_ref_discontinuous_cds.cpp)

add_executable(check_codon_extraction test_codon_extraction.cpp)

add_executable(check_write_codon_aa test_write_codon_aa.cpp)

set(ALL_TEST_EXECUTABLES check_primer_trim
        check_trim
        check_quality_trim
        check_consensus
        check_allele_depth
        check_consensus_threshold
        check_consensus_min_depth
        check_consensus_seq_id
        check_primer_bed
        check_clustering_consensus
        check_clustering_majority
        check_unit_consensus
        check_getmasked
        check_removereads
        check_variants
#        check_variants_2
        check_common_variants
        check_unpaired_trim
        check_primer_trim_edge_cases
        check_isize_trim
        check_interval_tree
        check_amplicon_search
        check_consensus_min_insert_threshold
        check_variant_caller
        check_site_aggregator
        check_site_amplicon_aggregator_stats
        check_site_state_aggregator_stats
        check_gmm_1d
        check_gaussian_sampler
        check_multinomial_sampler
        check_bootstrap_variants
        check_bd_threshold
        check_subset_sum
        check_parse_gff
        check_ref_discontinuous_cds
        check_codon_extraction
        check_write_codon_aa
)

# All shared code lives in ivar_core, which PUBLIC-propagates the htslib
# include directory and link dependencies to every test executable.
foreach(test_executable ${ALL_TEST_EXECUTABLES})
    target_link_libraries(${test_executable} PRIVATE ivar_core)
    add_test(NAME ${test_executable} COMMAND ${test_executable})
endforeach ()
