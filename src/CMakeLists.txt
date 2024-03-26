# Add library or executable target
file(GLOB armadillo_src
        "include/armadillo_bits/*.hpp"
)

add_executable(ivar
        ivar.cpp
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
        parse_gff.cpp
        ref_seq.cpp
        interval_tree.cpp
        saga.cpp
        gmm.cpp
        population_estimate.cpp
        call_consensus_clustering.cpp
        ${armadillo_src}
)

target_link_libraries(ivar z pthread)

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
    target_link_libraries(ivar ${HTSLIB})
endif()

if(HTSLIB_PREFIX)
    target_include_directories(ivar PRIVATE ${HTSLIB_PREFIX}/include)
else()
    message(FATAL_ERROR "htslib header files not found")
endif()

# Specify include directories if needed
#target_include_directories(ivar PUBLIC ${PROJECT_SOURCE_DIR}/include)