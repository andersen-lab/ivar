add_executable(check_primer_trim
        test_primer_trim.cpp
        ../src/trim_primer_quality.cpp
        ../src/primer_bed.cpp
        ../src/interval_tree.cpp
        ../src/allele_functions.cpp
        ../src/ref_seq.cpp
        ../src/parse_gff.cpp
)

add_executable(check_trim
        test_trim.cpp
        ../src/trim_primer_quality.cpp
        ../src/primer_bed.cpp
        ../src/interval_tree.cpp
        ../src/allele_functions.cpp
        ../src/ref_seq.cpp
        ../src/parse_gff.cpp
)

add_executable(check_quality_trim
        test_quality_trim.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_consensus
        test_call_consensus_from_plup.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
)

add_executable(check_allele_depth
        test_allele_depth.cpp ../src/allele_functions.cpp
)

add_executable(check_consensus_threshold
        test_consensus_threshold.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
)

add_executable(check_consensus_min_insert_threshold
        test_consensus_min_insert_threshold.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
)

add_executable(check_consensus_min_depth
        test_consensus_min_depth.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
)

add_executable(check_consensus_seq_id
        test_consensus_seq_id.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
)

add_executable(check_variants
        test_variants.cpp ../src/call_variants.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_common_variants
        test_common_variants.cpp ../src/get_common_variants.cpp
)

add_executable(check_primer_bed
        test_primer_bed.cpp ../src/primer_bed.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_clustering_consensus
        test_clustering_consensus.cpp ../src/gmm.cpp ../src/call_consensus_clustering.cpp ../src/estimate_error.cpp ../src/saga.cpp ../src/interval_tree.cpp ../src/primer_bed.cpp ../src/allele_functions.cpp ../src/trim_primer_quality.cpp
)

add_executable(check_getmasked
        test_getmasked.cpp ../src/get_masked_amplicons.cpp ../src/primer_bed.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_removereads
        test_removereads.cpp ../src/remove_reads_from_amplicon.cpp ../src/primer_bed.cpp ../src/trim_primer_quality.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_unpaired_trim
        test_unpaired_trim.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_primer_trim_edge_cases
        test_primer_trim_edge_cases.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_isize_trim
        test_isize_trim.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_interval_tree
        test_interval_tree.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_amplicon_search
        test_amplicon_search.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/allele_functions.cpp ../src/parse_gff.cpp ../src/ref_seq.cpp
)

add_executable(check_preprocessing
        test_preprocessing.cpp ../src/allele_functions.cpp ../src/saga.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp ../src/interval_tree.cpp ../src/ref_seq.cpp ../src/parse_gff.cpp
)

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
        check_getmasked
        check_removereads
        check_variants
        check_common_variants
        check_unpaired_trim
        check_primer_trim_edge_cases
        check_isize_trim
        check_interval_tree
        check_amplicon_search
        check_consensus_min_insert_threshold
        check_preprocessing
)

if(HTSLIB_PREFIX)
    find_library(HTSLIB
            NAMES hts
            PATHS ${HTSLIB_PREFIX}/lib
            NO_DEFAULT_PATH)
else()
    find_library(HTSLIB
            NAMES hts)
endif()

if(NOT HTSLIB)
    message(FATAL_ERROR "htslib not found")
else()
    foreach(test_executable ${ALL_TEST_EXECUTABLES})
        # Include directories for htslib
        target_include_directories(${test_executable} PRIVATE ${HTSLIB_PREFIX}/include)

        # Link the executable with htslib
        target_link_libraries(${test_executable} ${HTSLIB})

        add_test(NAME ${test_executable} COMMAND ${test_executable})
    endforeach ()
endif()
