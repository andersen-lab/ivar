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
endif()

if(NOT HTSLIB_PREFIX)
    message(FATAL_ERROR "htslib header files not found")
endif()

add_library(ivar_core STATIC
        call_consensus_pileup.cpp
        alignment.cpp
        suffix_tree.cpp
        trim_primer_quality.cpp
        remove_reads_from_amplicon.cpp
        call_variants.cpp
        primer_bed.cpp
        allele_functions.cpp
        get_masked_amplicons.cpp
        get_common_variants.cpp
        genomic_position.cpp
        variant_caller.cpp
        parse_gff.cpp
        ref_seq.cpp
        interval_tree.cpp
        saga.cpp
        site_state.cpp
        site_aggregator.cpp
        site_aggregator_stats.cpp
        site_state_aggregator_stats.cpp
        site_amplicon_aggregator_stats.cpp
        gmm_1d.cpp
        gmm.cpp
        call_consensus_clustering.cpp
        solve_clustering.cpp
        bootstrap_variants.cpp
        multinomial_sampler.cpp
        gaussian_sampler.cpp
)

target_link_libraries(ivar_core PUBLIC z pthread ${HTSLIB})
target_include_directories(ivar_core PUBLIC ${HTSLIB_PREFIX}/include)

add_executable(ivar ivar.cpp)
target_link_libraries(ivar PRIVATE ivar_core)
