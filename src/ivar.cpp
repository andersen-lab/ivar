/*! \file */

#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "call_consensus_pileup.h"
#include "call_variants.h"
#include "get_common_variants.h"
#include "get_masked_amplicons.h"
#include "remove_reads_from_amplicon.h"
#include "suffix_tree.h"
#include "trim_primer_quality.h"
#include "gmm.h"
#include "call_consensus_clustering.h"
#include "saga.h"

const std::string VERSION = "1.4.2";

struct args_t {
  std::string bam;               // -i
  std::string bed;               // -b
  std::string text;              // -t
  std::string seq_id;            // -i for consensus
  std::string prefix;            // -p
  std::string ref;               // -r
  std::string region;            // -R
  uint8_t min_qual;              // -q
  uint8_t sliding_window;        // -s
  double min_threshold;          // -t
  double min_insert_threshold;   // -c
  int min_length;                // -m
  std::string f1;                // -1
  std::string f2;                // -2
  std::string adp_path;          // -a
  uint32_t min_depth;             // -m
  char gap;                      // -n
  bool keep_min_coverage;        // -k
  std::string primer_pair_file;  // -f
  int32_t primer_offset;         // -x
  std::string file_list;         // -f
  bool write_no_primers_flag;    // -e
  std::string gff;               // -g
  bool keep_for_reanalysis;      // -k
  //contam params
  std::string variants;          // -s
  float dcov_1;                  // -d
  float dcov_2;                  // -D
} g_args;

void print_usage() {
  std::cout
      << "Usage:	ivar [command "
         "<trim|variants|filtervariants|consensus|getmasked|removereads|contam|saga|"
         "version|help>]\n"
         "\n"
         "        Command       Description\n"
         "           trim       Trim reads in aligned BAM file\n"
         "         contam       Estimate the number of populations in a sample\n"
         "           saga       Call consensus with an automated threshold\n"
         "       variants       Call variants from aligned BAM file\n"
         " filtervariants       Filter variants across replicates or samples\n"
         "      consensus       Call consensus from aligned BAM file\n"
         "      getmasked       Detect primer mismatches and get primer "
         "indices for the amplicon to be masked\n"
         "    removereads       Remove reads from trimmed BAM file\n"
         "        version       Show version information\n"
         "\n"
         "To view detailed usage for each command type `ivar <command>` \n";
}

void print_contam_usage() {
  std::cout
      << "Usage: ivar contam \n\n"
         "Input Options    Description\n"
         "           -i    BAM file, with aligned reads, to "
         "trim primers and quality. If not specified will use standard in\n"
         "           -b    BED file with primer sequences and positions. If no "
         "BED file is specified, only quality trimming will be done.\n"
         "           -f    [EXPERIMENTAL] Primer pair information file "
         "containing left and right primer names for the same amplicon "
         "separated by a tab\n"
         "                 If provided, reads that do not fall within atleat "
         "one amplicon will be ignored prior to primer trimming.\n"
         "           -x    Primer position offset (Default: 0). Reads that "
         "occur at the specified offset positions relative to primer positions "
         "will also be trimmed.\n"
         "           -m    Minimum length of read to retain after trimming "
         "(Default: 50% average length of the first 1000 reads)\n"
         "           -q    Minimum quality threshold for sliding window to "
         "pass (Default: 20)\n"
         "           -s    Width of sliding window (Default: 4)\n"
         "           -e    Include reads with no primers. By default, reads "
         "with no primers are excluded\n"
         "           -k    Keep reads to allow for reanalysis: keep reads "
         "which would be dropped by\n"
         "                 alignment length filter or primer requirements, but "
         "mark them QCFAIL\n\n"
         "Output Options   Description\n"
         "           -p    Prefix for the output BAM file. If none is specified output will go to std out\n";
}

void print_trim_usage() {
  std::cout
      << "Usage: ivar trim -i [<input.bam>] -b <primers.bed> [-p <prefix>] [-m "
         "<min-length>] [-q <min-quality>] [-s <sliding-window-width>]\n\n"
         "Input Options    Description\n"
         "           -i    BAM file, with aligned reads, to "
         "trim primers and quality. If not specified will use standard in\n"
         "           -b    BED file with primer sequences and positions. If no "
         "BED file is specified, only quality trimming will be done.\n"
         "           -f    [EXPERIMENTAL] Primer pair information file "
         "containing left and right primer names for the same amplicon "
         "separated by a tab\n"
         "                 If provided, reads that do not fall within atleat "
         "one amplicon will be ignored prior to primer trimming.\n"
         "           -x    Primer position offset (Default: 0). Reads that "
         "occur at the specified offset positions relative to primer positions "
         "will also be trimmed.\n"
         "           -m    Minimum length of read to retain after trimming "
         "(Default: 50% average length of the first 1000 reads)\n"
         "           -q    Minimum quality threshold for sliding window to "
         "pass (Default: 20)\n"
         "           -s    Width of sliding window (Default: 4)\n"
         "           -e    Include reads with no primers. By default, reads "
         "with no primers are excluded\n"
         "           -k    Keep reads to allow for reanalysis: keep reads "
         "which would be dropped by\n"
         "                 alignment length filter or primer requirements, but "
         "mark them QCFAIL\n\n"
         "Output Options   Description\n"
         "           -p    Prefix for the output BAM file. If none is specified output will go to std out\n";
}

void print_variants_usage() {
  std::cout
      << "Usage: samtools mpileup -aa -A -d 0 -B -Q 0 --reference "
         "[<reference-fasta] <input.bam> | ivar variants -p <prefix> [-q "
         "<min-quality>] [-t <min-frequency-threshold>] [-m <minimum depth>] "
         "[-r <reference-fasta>] [-g GFF file]\n\n"
         "Note : samtools mpileup output must be piped into ivar variants\n\n"
         "Input Options    Description\n"
         "           -q    Minimum quality score threshold to count base "
         "(Default: 20)\n"
         "           -t    Minimum frequency threshold(0 - 1) to call variants "
         "(Default: 0.03)\n"
         "           -m    Minimum read depth to call variants (Default: 0)\n"
         "           -r    Reference file used for alignment. This is used to "
         "translate the nucleotide sequences and identify intra host single "
         "nucleotide variants\n"
         "           -g    A GFF file in the GFF3 format can be supplied to "
         "specify coordinates of open reading frames (ORFs). In absence of GFF "
         "file, amino acid translation will not be done.\n\n"
         "Output Options   Description\n"
         "           -p    (Required) Prefix for the output tsv variant "
         "file\n\n";
}

void print_filtervariants_usage() {
  std::cout
      << "Usage: ivar filtervariants -p <prefix> replicate-one.tsv "
         "replicate-two.tsv ... OR ivar filtervariants -p <prefix> -f <text "
         "file with one variant file per line> \n"
         "Input: Variant tsv files for each replicate/sample\n\n"
         "Input Options    Description\n"
         "           -t    Minimum fraction of files required to contain the "
         "same variant. Specify value within [0,1]. (Default: 1)\n"
         "           -f    A text file with one variant file per line.\n\n"
         "Output Options   Description\n"
         "           -p    (Required) Prefix for the output filtered tsv "
         "file\n";
}

void print_consensus_usage() {
  std::cout
      << "Usage: samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar "
         "consensus -p <prefix> \n\n"
         "Note : samtools mpileup output must be piped into `ivar "
         "consensus`\n\n"
         "Input Options    Description\n"
         "           -q    Minimum quality score threshold to count base "
         "(Default: 20)\n"
         "           -t    Minimum frequency threshold(0 - 1) to call "
         "consensus. (Default: 0)\n"
         "                 Frequently used thresholds | Description\n"
         "                 ---------------------------|------------\n"
         "                                          0 | Majority or most "
         "common base\n"
         "                                        0.2 | Bases that make up "
         "atleast 20% of the depth at a position\n"
         "                                        0.5 | Strict or bases that "
         "make up atleast 50% of the depth at a position\n"
         "                                        0.9 | Strict or bases that "
         "make up atleast 90% of the depth at a position\n"
         "                                          1 | Identical or bases "
         "that make up 100% of the depth at a position. Will have highest "
         "ambiguities\n"
         "           -c    Minimum insertion frequency threshold(0 - 1) to "
         "call consensus. (Default: 0.8)\n"
         "                 Frequently used thresholds | Description\n"
         "                 ---------------------------|------------\n"
         "                                          0 | Allow insertion if it "
         "appears even once\n"
         "                                        0.2 | Insertions with at "
         "least 20% of the depth at a position\n"
         "                                        0.5 | Insertion with at "
         "least 50% of the depth at a position\n"
         "                                        0.9 | Insertions with at "
         "least 90% of the depth at a position\n"
         "                                          1 | Insertion with 100% of "
         "the depth at a position. Will have highest ambiguities\n"
         "           -m    Minimum depth to call consensus(Default: 10)\n"
         "           -k    If '-k' flag is added, regions with depth less than "
         "minimum depth will not be added to the consensus sequence. Using "
         "'-k' will override any option specified using -n \n"
         "           -n    (N/-) Character to print in regions with less than "
         "minimum coverage(Default: N)\n\n"
         "Output Options   Description\n"
         "           -p    (Required) Prefix for the output fasta file and "
         "quality file\n"
         "           -i    (Optional) Name of fasta header. By default, the "
         "prefix is used to create the fasta header in the following format, "
         "Consensus_<prefix>_threshold_<frequency-threshold>_quality_<minimum-"
         "quality>_<min-insert-threshold>\n";
}

void print_removereads_usage() {
  std::cout
      << "Usage: ivar removereads -i <input.trimmed.bam> -p <prefix> -t "
         "<text-file-with-primer-indices> -b <primers.bed> \n"
         "Note: This step is used only for amplicon-based sequencing.\n\n"
         "Input Options    Description\n"
         "           -i    (Required) Input BAM file  trimmed with ‘ivar "
         "trim’. Must be sorted which can be done using `samtools sort`.\n"
         "           -t    (Required) Text file with primer indices separated "
         "by spaces. This is the output of `getmasked` command.\n"
         "           -b    (Required) BED file with primer sequences and "
         "positions.\n\n"
         "Output Options   Description\n"
         "           -p    (Required) Prefix for the output filtered BAM "
         "file\n";
}

void print_getmasked_usage() {
  std::cout << "Usage: ivar getmasked -i <input-filtered.tsv> -b <primers.bed> "
               "-f <primer_pairs.tsv> -p <prefix>\n"
               "Note: This step is used only for amplicon-based sequencing.\n\n"
               "Input Options    Description\n"
               "           -i    (Required) Input filtered variants tsv "
               "generated from `ivar filtervariants`\n"
               "           -b    (Required) BED file with primer sequences and "
               "positions\n"
               "           -f    (Required) Primer pair information file "
               "containing left and right primer names for the same amplicon "
               "separated by a tab\n"
               "Output Options   Description\n"
               "           -p    (Required) Prefix for the output text file\n";
}

void print_trimadapter_usage() {
  std::cout << "NOTE: EXPERIMENTAL FEATURE\n"
               "Usage: ivar trimadapter [-f1 <input-fastq>] [-f2 "
               "<input-fastq-2>] [-p prefix] [-a <adapter-fasta-file>]\n\n"
               "Input Options    Description\n"
               "           -1    (Required) Input fastq file\n"
               "           -2    Input fastq file 2 (for pair ended reads)\n"
               "           -a    (Required) Adapter Fasta File\n\n"
               "Output Options   Description\n"
               "           -p    (Required) Prefix of output fastq files\n";
}

void print_version_info() {
  std::cout << "iVar version " << VERSION << std::endl
            << "\nPlease raise issues and bug reports at "
               "https://github.com/andersen-lab/ivar/\n\n";
}

static const char *trim_opt_str = "i:b:f:x:p:m:q:s:ekh?";
static const char *variants_opt_str = "p:t:q:m:r:g:h?";
static const char *consensus_opt_str = "i:p:q:t:c:m:n:kh?";
static const char *removereads_opt_str = "i:p:t:b:h?";
static const char *filtervariants_opt_str = "p:t:f:h?";
static const char *getmasked_opt_str = "i:b:f:p:h?";
static const char *trimadapter_opt_str = "1:2:p:a:h?";
static const char *contam_opt_str = "i:b:f:x:p:m:q:s:d:D:e:r:kh?";

std::string get_filename_without_extension(std::string f, std::string ext) {
  if (ext.length() > f.length())  // If extension longer than filename
    return f;
  if (f.substr(f.length() - ext.length(), ext.length()).compare(ext) == 0) {
    return f.substr(0, f.length() - ext.length());
  }
  return f;
}

/*!
 Main Function

 This is where the command line arguments into iVar are parsed.
 iVar first parses the first argument and depending on the command it either
 returns the version or uses GNU getopt to parse the remaining arguments.
 */

int main(int argc, char *argv[]) {
  if (argc == 1) {
    print_usage();
    return -1;
  }
  std::stringstream cl_cmd;
  cl_cmd << "@PG\tID:ivar-" << argv[1] << "\tPN:ivar\tVN:" << VERSION
         << "\tCL:" << argv[0] << " ";
  for (int i = 1; i < argc; ++i) {
    cl_cmd << argv[i];
    if (i != argc - 1) cl_cmd << " ";
  }
  cl_cmd << "\n\0";
  std::string cmd(argv[1]);
  if (cmd.compare("-v") == 0) {
    print_version_info();
    return 0;
  }
  int opt = 0, res = 0;
  // Sift arg by 1 for getopt
  argv[1] = argv[0];
  argv++;
  argc--;

  //ivar contam
  if (cmd.compare("contam") == 0) {
    g_args.bed = "";
    g_args.primer_pair_file = "";
    g_args.primer_offset = 0;
    opt = getopt(argc, argv, contam_opt_str);
    while (opt != -1) {
      switch (opt) {
        case 'i':
          g_args.bam = optarg;
          break;
        case 'b':
          g_args.bed = optarg;
          break;
        case 'f':
          g_args.primer_pair_file = optarg;
          break;
        case 'x':
          g_args.primer_offset = std::stoi(optarg);
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 's':
          g_args.variants = optarg;
          break;
        case 'd':
          g_args.dcov_1 = std::stof(optarg);
          break;
        case 'D':
          g_args.dcov_2 = std::stof(optarg);
          break;
        case 'h':
        case '?':
          print_trim_usage();
          return -1;
          break;
      }
      opt = getopt(argc, argv, contam_opt_str);
    }
    if (!g_args.variants.empty() && !g_args.prefix.empty()) {
      std::vector<variant> variants = gmm_model(g_args.variants, g_args.prefix);
      cluster_consensus(variants, g_args.prefix);
    }
    res = 0; 
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".bam");
    //cluster_consensus();
  }

  //ivar saga
  else if (cmd.compare("saga") == 0) {
    g_args.bed = "";
    g_args.primer_pair_file = "";
    g_args.primer_offset = 0;
    opt = getopt(argc, argv, trim_opt_str);
    while (opt != -1) {
      switch (opt) {
        case 'i':
          g_args.bam = optarg;
          break;
        case 'b':
          g_args.bed = optarg;
          break;
        case 'f':
          g_args.primer_pair_file = optarg;
          break;
        case 'x':
          g_args.primer_offset = std::stoi(optarg);
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 'm':
          g_args.min_length = std::stoi(optarg);
          break;
        case 'h':
        case '?':
          print_trim_usage();
          return -1;
          break;
      }
      opt = getopt(argc, argv, trim_opt_str);
    }
    if (g_args.bam.empty() && isatty(STDIN_FILENO)) {
      std::cout << "Please supply a BAM file using -i or supply the input file "
                   "through standard input"
                << std::endl
                << std::endl;
      print_trim_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".bam");
    res = preprocess_reads(g_args.bam, g_args.bed, g_args.prefix,
                               cl_cmd.str(),
                               g_args.primer_pair_file, g_args.primer_offset);
    //std::vector<uint32_t> populations_iterate;
    //for(uint32_t i= 2; i < 6; i++){
    //  populations_iterate.push_back(i);
    //}
    //res = gmm_model(g_args.prefix, populations_iterate, g_args.prefix);
    //res = cluster_consensus(
  }

  // ivar trim
  else if (cmd.compare("trim") == 0) {
    g_args.min_qual = 20;
    g_args.sliding_window = 4;
    g_args.min_length = -1;
    g_args.write_no_primers_flag = false;
    g_args.keep_for_reanalysis = false;
    g_args.bed = "";
    g_args.primer_pair_file = "";
    g_args.primer_offset = 0;
    opt = getopt(argc, argv, trim_opt_str);
    while (opt != -1) {
      switch (opt) {
        case 'i':
          g_args.bam = optarg;
          break;
        case 'b':
          g_args.bed = optarg;
          break;
        case 'f':
          g_args.primer_pair_file = optarg;
          break;
        case 'x':
          g_args.primer_offset = std::stoi(optarg);
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 'm':
          g_args.min_length = std::stoi(optarg);
          break;
        case 'q':
          g_args.min_qual = std::stoi(optarg);
          break;
        case 's':
          g_args.sliding_window = std::stoi(optarg);
          break;
        case 'e':
          g_args.write_no_primers_flag = true;
          break;
        case 'k':
          g_args.keep_for_reanalysis = true;
          break;
        case 'h':
        case '?':
          print_trim_usage();
          return -1;
          break;
      }
      opt = getopt(argc, argv, trim_opt_str);
    }
    if (g_args.bam.empty() && isatty(STDIN_FILENO)) {
      std::cout << "Please supply a BAM file using -i or supply the input file "
                   "through standard input"
                << std::endl
                << std::endl;
      print_trim_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".bam");
    res = trim_bam_qual_primer(g_args.bam, g_args.bed, g_args.prefix,
                               g_args.min_qual, g_args.sliding_window,
                               cl_cmd.str(), g_args.write_no_primers_flag,
                               g_args.keep_for_reanalysis, g_args.min_length,
                               g_args.primer_pair_file, g_args.primer_offset);
  }
  // ivar variants
  else if (cmd.compare("variants") == 0) {
    g_args.min_qual = 20;
    g_args.min_threshold = 0.03;
    g_args.min_depth = 0;
    g_args.ref = "";
    g_args.gff = "";
    opt = getopt(argc, argv, variants_opt_str);
    while (opt != -1) {
      switch (opt) {
        case 'p':
          g_args.prefix = optarg;
          break;
        case 't':
          g_args.min_threshold = atof(optarg);
          break;
        case 'q':
          g_args.min_qual = std::stoi(optarg);
          break;
        case 'm':
          g_args.min_depth = std::stoi(optarg);
          break;
        case 'r':
          g_args.ref = optarg;
          break;
        case 'g':
          g_args.gff = optarg;
          break;
        case 'h':
        case '?':
          print_variants_usage();
          return 0;
          break;
      }
      opt = getopt(argc, argv, variants_opt_str);
    }
    if (g_args.prefix.empty()) {
      print_variants_usage();
      return -1;
    }
    if (g_args.gff.empty())
      std::cout << "A GFF file containing the open reading frames (ORFs) has "
                   "not been provided. Amino acid translation will not be done."
                << std::endl;
    if (g_args.ref.empty())
      std::cout << "A reference sequence has not been supplied. Amino acid "
                   "translation will not be done."
                << std::endl;
    if (!g_args.gff.empty() &&
        g_args.ref.empty()) {  // GFF specified but no reference then exit
      std::cout << "Please specify reference (using -r) based on which the GFF "
                   "file was computed."
                << std::endl;
      print_variants_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".tsv");
    g_args.min_threshold =
        (g_args.min_threshold < 0 || g_args.min_threshold > 1)
            ? 0.03
            : g_args.min_threshold;
    if (isatty(STDIN_FILENO)) {
      std::cout << "Please pipe mpileup into `ivar variants` command.\n\n";
      print_variants_usage();
      return -1;
    }
    res = call_variants_from_plup(std::cin, g_args.prefix, g_args.min_qual,
                                  g_args.min_threshold, g_args.min_depth,
                                  g_args.ref, g_args.gff);
  }
  // ivar consensus
  else if (cmd.compare("consensus") == 0) {
    opt = getopt(argc, argv, consensus_opt_str);
    g_args.seq_id = "";
    g_args.min_threshold = 0;
    g_args.min_depth = 10;
    g_args.gap = 'N';
    g_args.min_qual = 20;
    g_args.keep_min_coverage = true;
    g_args.min_insert_threshold = 0.8;
    while (opt != -1) {
      switch (opt) {
        case 't':
          g_args.min_threshold = atof(optarg);
          break;
        case 'c':
          g_args.min_insert_threshold = atof(optarg);
          break;
        case 'i':
          g_args.seq_id = optarg;
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 'm':
          g_args.min_depth = std::stoi(optarg);
          break;
        case 'n':
          g_args.gap = optarg[0];
          break;
        case 'q':
          g_args.min_qual = std::stoi(optarg);
          break;
        case 'k':
          g_args.keep_min_coverage = false;
        case 'g':
          break;
        case 'h':
        case '?':
          print_consensus_usage();
          return 0;
          break;
      }
      opt = getopt(argc, argv, consensus_opt_str);
    }
    if (g_args.prefix.empty()) {
      print_consensus_usage();
      return -1;
    }
    if (isatty(STDIN_FILENO)) {
      std::cout << "Please pipe mpileup into `ivar consensus` command.\n\n";
      print_consensus_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".fa");
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".fasta");
    g_args.gap = (g_args.gap != 'N' && g_args.gap != '-')
                     ? 'N'
                     : g_args.gap;  // Accept only N or -
    std::cout << "Minimum Quality: " << (uint16_t)g_args.min_qual << std::endl;
    std::cout << "Threshold: " << g_args.min_threshold << std::endl;
    std::cout << "Minimum depth: " << (unsigned)g_args.min_depth << std::endl;
    std::cout << "Minimum Insert Threshold: " << g_args.min_insert_threshold
              << std::endl;
    if (!g_args.keep_min_coverage)
      std::cout << "Regions with depth less than minimum depth will not added "
                   "to consensus"
                << std::endl;
    else
      std::cout << "Regions with depth less than minimum depth covered by: "
                << g_args.gap << std::endl;
    res = call_consensus_from_plup(
        std::cin, g_args.seq_id, g_args.prefix, g_args.min_qual,
        g_args.min_threshold, g_args.min_depth, g_args.gap,
        g_args.keep_min_coverage, g_args.min_insert_threshold);
  } else if (cmd.compare("removereads") == 0) {
    opt = getopt(argc, argv, removereads_opt_str);
    while (opt != -1) {
      switch (opt) {
        case 'i':
          g_args.bam = optarg;
          break;
        case 't':
          g_args.text = optarg;
          break;
        case 'b':
          g_args.bed = optarg;
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 'h':
        case '?':
          print_removereads_usage();
          return 0;
          break;
      }
      opt = getopt(argc, argv, removereads_opt_str);
    }
    if (g_args.bam.empty() || g_args.prefix.empty() || g_args.bed.empty() ||
        g_args.text.empty()) {
      print_removereads_usage();
      return -1;
    }
    std::string s;
    std::vector<std::string> amp;
    std::ifstream fin(g_args.text.c_str());
    while (getline(fin, s, '\t')) {
      amp.push_back(s);
    }
    fin.close();
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".bam");
    res = rmv_reads_from_amplicon(g_args.bam, g_args.region, g_args.prefix, amp,
                                  g_args.bed, cl_cmd.str());

  } else if (cmd.compare("filtervariants") == 0) {
    opt = getopt(argc, argv, filtervariants_opt_str);
    g_args.min_threshold = 1;
    while (opt != -1) {
      switch (opt) {
        case 'p':
          g_args.prefix = optarg;
          break;
        case 't':
          g_args.min_threshold = atof(optarg);
          break;
        case 'f':
          g_args.file_list = optarg;
          break;
        case 'h':
        case '?':
          print_filtervariants_usage();
          return 0;
          break;
      }
      opt = getopt(argc, argv, filtervariants_opt_str);
    }
    if (g_args.min_threshold < 0 || g_args.min_threshold > 1) {
      print_filtervariants_usage();
      return -1;
    }
    if (optind >= argc && g_args.file_list.empty()) {
      print_filtervariants_usage();
      return -1;
    }
    if (g_args.prefix.empty()) {
      print_filtervariants_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".tsv");
    // Read files from list
    char **files = new char *[100];
    int nfiles = 100, ctr = 0;
    std::string line;
    if (!g_args.file_list.empty()) {  // File list supplied
      std::ifstream file_fin = std::ifstream(g_args.file_list);
      while (std::getline(file_fin, line)) {
        files[ctr] = strdup(line.c_str());
        if (ctr == nfiles - 1) {
          nfiles += 100;
          *files = (char *)realloc(*files, nfiles * (sizeof(char *)));
        }
        ctr++;
      }
      file_fin.close();
      nfiles = (nfiles > ctr) ? ctr : nfiles;
      res = common_variants(g_args.prefix, g_args.min_threshold, files, nfiles);
      // Free files, nfiles
      for (int i = 0; i < nfiles; ++i) {
        free(files[i]);
      }
      delete[] files;
    } else {
      res = common_variants(g_args.prefix, g_args.min_threshold, argv + optind,
                            argc - optind);
    }
  } else if (cmd.compare("getmasked") == 0) {
    opt = getopt(argc, argv, getmasked_opt_str);
    while (opt != -1) {
      switch (opt) {
        case 'i':
          g_args.bam = optarg;
          break;
        case 'b':
          g_args.bed = optarg;
          break;
        case 'f':
          g_args.primer_pair_file = optarg;
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 'h':
        case '?':
          print_getmasked_usage();
          return 0;
          break;
      }
      opt = getopt(argc, argv, getmasked_opt_str);
    }
    if (g_args.bed.empty() || g_args.bam.empty() || g_args.prefix.empty() ||
        g_args.primer_pair_file.empty()) {
      print_getmasked_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix, ".txt");
    res = get_primers_with_mismatches(g_args.bed, g_args.bam, g_args.prefix,
                                      g_args.primer_pair_file);
  } else if (cmd.compare("trimadapter") == 0) {
    opt = getopt(argc, argv, trimadapter_opt_str);
    while (opt != -1) {
      switch (opt) {
        case '1':
          g_args.f1 = optarg;
          break;
        case '2':
          g_args.f2 = optarg;
          break;
        case 'p':
          g_args.prefix = optarg;
          break;
        case 'a':
          g_args.adp_path = optarg;
          break;
        case 'h':
        case '?':
          print_trimadapter_usage();
          return 0;
          break;
      }
      opt = getopt(argc, argv, trimadapter_opt_str);
    }
    if (g_args.f1.empty() || g_args.prefix.empty() || g_args.adp_path.empty()) {
      print_trimadapter_usage();
      return -1;
    }
    res = trim_adapter(g_args.f1, g_args.f2, g_args.adp_path, g_args.prefix);
  } else if (cmd.compare("version") == 0) {
    print_version_info();
  } else {
    std::cout << "Unknown command: \"" << cmd << "\"" << std::endl << std::endl;
    print_usage();
  }
  return res;
}
