Manual {#manualpage}
============

[TOC]

Available Commands
====
| Command | Description |
|:--------|:------------|
| trim | Trim reads in aligned BAM |
| variants | Call variants from aligned BAM file |
| filtervariants | Filter variants across replicates or multiple samples aligned using the same reference |
| consensus | Call consensus from aligned BAM file |
| getmasked | Detect primer mismatches and get primer indices for the amplicon to be masked |
| removereads | Remove reads from trimmed BAM file |
| version | Show version information |
| trimadapter | (EXPERIMENTAL) Trim adapter sequences from reads |

To view detailed usage for each command type `ivar <command>`
Note : Commands maked (EXPERIMENTAL) are still under active development.

Description of commands
====

Trim primer sequences with iVar
----

iVar uses primer positions supplied in a BED file to soft clip primer sequences from an aligned and sorted BAM file. Following this, the reads are trimmed based on a quality threshold(Default: 20). To do the quality trimming, iVar uses a sliding window approach(Default: 4). The windows slides from the 5' end to the 3' end and if at any point the average base quality in the window falls below the threshold, the remaining read is soft clipped. If after trimming, the length of the read is greater than the minimum length specified(Default: 30), the read is written to the new trimmed BAM file.

Please note that the strand is taken into account while doing the trimming so forward primers are trimmed only from forward strand and reverse primers are trimmed from reverse strand.

To sort and index an aligned BAM file (OPTIONAL, if index is not present iVar will create one), the following command can be used,

```
# Input file - test.bam
samtools sort -o test.sorted.bam test.bam && samtools index test.sorted.bam.
```

**Note**: All the trimming in iVar is done by soft-clipping reads in an aligned BAM file. This information is lost if reads are extracted in fastq or fasta format from the trimmed BAM file.

Command:
```
ivar trim

Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]

Input Options    Description
           -i    (Required) Sorted bam file, with aligned reads, to trim primers and quality
           -b    (Required) BED file with primer sequences and positions
           -f    Primer pair information file containing left and right primer names for the same amplicon separated by a tab
                 If provided, reads will be filtered based on their overlap with amplicons prior to trimming
           -m    Minimum length of read to retain after trimming (Default: 30)
           -q    Minimum quality threshold for sliding window to pass (Default: 20)
           -s    Width of sliding window (Default: 4)
           -e    Include reads with no primers. By default, reads with no primers are excluded
           -k    Keep reads to allow for reanalysis: keep reads which would be dropped by
                 alignment length filter or primer requirements, but mark them QCFAIL

Output Options   Description
           -p    (Required) Prefix for the output BAM file
```

Example Usage:
```
ivar trim -b test_primers.bed -p test.trimmed -i test.bam -q 15 -m 50 -s 4
```

The command above will produce a trimmed BAM file test.trimmed.bam after trimming the aligned reads in test.bam using the primer positions specified in test_primers.bed and a minimum quality threshold of **15**, minimum read length of **50** and a sliding window of **4**.

Example BED file

```
Puerto	28	52	400_1_out_L	60	+
Puerto	482	504	400_1_out_R	60	-
Puerto	359	381	400_2_out_L	60	+
Puerto	796	818	400_2_out_R	60	-
Puerto	658	680	400_3_out_L*	60	+
Puerto	1054	1076	400_3_out_R*	60	-
.
.
.
.
```

Call variants with iVar
----

iVar uses the output of the `samtools mpileup` command to call variants - single nucleotide variants(SNVs) and indels. In order to call variants correctly, the reference file used for alignment must be passed to iVar using the `-r` flag. The output of `samtools pileup` is piped into `ivar variants` to generate a .tsv file with the variants. There are two parameters that can be set for variant calling using iVar - minimum quality(Default: 20) and minimum frequency(Default: 0.03). Minimum quality is the minimum quality for a base to be counted towards the ungapped depth to canculate iSNV frequency at a given position. For insertions, the quality metric is discarded and the mpileup depth is used directly. Minimum frequency is the minimum frequency required for a SNV or indel to be reported. 

#### Amino acid translation of iSNVs

iVar can identify codons and translate variants into amino acids using a GFF file in the [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) format containing the required coding regions (CDS). In absence of a GFF file, iVar will not perform the translation and "NA" will be added to the output file in place of the reference and alternate codons and amino acids. The GFF file in the GFF3 format can be downloaded via ftp from NCBI RefSeq/Genbank. They are usually the files with the extension ".gff.gz". For example, the GFF file for Zaire Ebolavirus can be found [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Zaire_ebolavirus/all_assembly_versions/GCF_000848505.1_ViralProj14703). More details on GFF3 files hosted by NCBI can be found in their ftp [FAQs](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/).

#### Account for RNA editing through polymerase slippage

Some RNA viruses such as Ebola virus, might have polymerase slippage causing the insertion of a couple of nucleotides. More details can be found [here](https://viralzone.expasy.org/857?outline=all_by_protein). iVar can account for this editing and identify the correct open reading frames. The user will have to specify two additional parameters, **EditPosition**: Position at which edit occurs and **EditSequence**: The sequence tht is inserted at the given positon, in the "attributes" column of the GFF file to account for this. A test example is given below,

```
test	Genbank	CDS	2	292	.	+	.	ID=id-testedit1;Note=PinkFloyd;EditPosition=100;EditSequence=A
test	Genbank	CDS	2	292	.	+	.	ID=id-testedit2;Note=AnotherBrickInTheWall;EditPosition=102;EditSequence=AA
```

If a certain base is present in multiple CDSs, iVar will add a new row for each CDS frame and distinguish the rows by adding the ID (specified in attributes of GFF) of the GFF feature used for the translation. This is shown for position 42 in the example output below. There are two rows with two different GFF features: id-test3 and id-test4. 

Command:
```
Usage: samtools mpileup -aa -A -d 0 -B -Q 0 --reference [<reference-fasta] <input.bam> | ivar variants -p <prefix> [-q <min-quality>] [-t <min-frequency-threshold>] [-m <minimum depth>] [-r <reference-fasta>] [-g GFF file]

Note : samtools mpileup output must be piped into ivar variants

Input Options    Description
           -q    Minimum quality score threshold to count base (Default: 20)
           -t    Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
           -m    Minimum read depth to call variants (Default: 0)
           -r    Reference file used for alignment. This is used to translate the nucleotide sequences and identify intra host single nucleotide variants
           -g    A GFF file in the GFF3 format can be supplied to specify coordinates of open reading frames (ORFs). In absence of GFF file, amino acid translation will not be done.

Output Options   Description
           -p    (Required) Prefix for the output tsv variant file
```

Example Usage:
```
samtools mpileup -aa -A -d 600000 -B -Q 0 test.trimmed.bam | ivar variants -p test -q 20 -t 0.03 -r test_reference.fa -g test.gff
```

The command above will generate a test.tsv file.

Example of output .tsv file.

```
REGION	POS	REF	ALT	REF_DP	REF_RV	REF_QUAL	ALT_DP	ALT_RV	ALT_QUAL	ALT_FREQ	TOTAL_DP	PVAL	PASS	GFF_FEATURE	REF_CODON	REF_AA	ALT_CODON	ALT_AA
test	42	G	T	0	0	0	1	0	49	1	1	1	FALSE	id-test3	AGG	R	ATG	M
test	42	G	T	0	0	0	1	0	49	1	1	1	FALSE	id-test4	CAG	Q	CAT	H
test	320	A	T	1	1	35	1	1	46	0.5	2	0.666667	FALSE	NA	NA	NA	NA	NA
test	365	A	T	0	0	0	1	1	27	1	1	1	FALSE	NA	NA	NA	NA	NA
```

Description

| Field     | Description                               |
|:----------|:------------------------------------------|
| REGION    | Region from BAM file                      |
| POS       | Position on reference sequence            |
| REF       | Reference base                            |
| ALT       | Alternate Base                            |
| REF\_DP   | Ungapped depth of reference base                   |
| REF\_RV   | Ungapped depth of reference base on reverse reads  |
| REF\_QUAL | Mean quality of reference base            |
| ALT\_DP   | Ungapped depth of alternate base.                   |
| ALT\_RV   | Ungapped deapth of alternate base on reverse reads |
| ALT\_QUAL | Mean quality of alternate base            |
| ALT\_FREQ | Frequency of alternate base               |
| TOTAL\_DP | Total depth at position                   |
| PVAL      | p-value of fisher's exact test            |
| PASS      | Result of p-value <= 0.05                 |
| GFF\_FEATURE | ID of the GFF feature used for the translation |
| REF\_CODON | Codong using the reference base |
| REF\_AA | Amino acid translated from reference codon |
| ALT\_CODON | Codon using the alternate base |
| ALT\_AA | Amino acid translated from the alternate codon |

**Note**: Please use the -B options with `samtools mpileup` to call variants and generate consensus. When a reference sequence is supplied, the quality of the reference base is reduced to 0 (ASCII: !) in the mpileup output. Disabling BAQ with -B seems to fix this. This was tested in samtools 1.7 and 1.8.

Filter variants across replicates with iVar
----

iVar can be used to get an intersection of variants(in .tsv files) called from any number of replicates or from different samples using the same reference sequence. This intersection will filter out any iSNVs that do not occur in a minimum fraction of the files supplied. This parameter can be changed using the `-t` flag which range from 0 to 1 (default). Fields that are different across replicates(fields apart from REGION, POS, REF, ALT, REF\_CODON, REF\_AA, ALT\_CODON, ALT\_AA) will have the filename added as a suffix. If there are a large number of files to be filtered, the `-f` flag can be used to supply a text file with one sample/replicate variant .tsv file per line.

Command:
```
Usage: ivar filtervariants -p <prefix> replicate-one.tsv replicate-two.tsv ... OR ivar filtervariants -p <prefix> -f <text file with one variant file per line> 
Input: Variant tsv files for each replicate/sample

Input Options    Description
           -t    Minimum fraction of files required to contain the same variant. Specify value within [0,1]. (Default: 1)
           -f    A text file with one variant file per line.

Output Options   Description
           -p    (Required) Prefix for the output filtered tsv file
```

Example Usage:
The command below only retains those variants that are found in atleast 50% of the fiels supplied
```
ivar filtervariants -t 0.5 -p test.filtered test.1.tsv test.2.tsv test.3.tsv
```

The three replicates can also be supplied using a text file as shown below

```
ivar filtervariants -t 0.5 -p test.filtered -f filter_files.txt
```

filter_files.txt
```
./path/to/test.1.tsv
./path/to/test.2.tsv
./path/to/test.3.tsv
```

The command above will prodoce an output .tsv file test.filtered.tsv.

Example output of filtered .tsv file from three files test_rep1.tsv and test_rep2.tsv

```
REGION	POS	REF	ALT	GFF_FEATURE	REF_CODON	REF_AA	ALT_CODON	ALT_AA	REF_DP_test.1.tsv	REF_RV_test.1.tsv	REF_QUAL_test.1.tsv	ALT_DP_test.1.tsv	ALT_RV_test.1.tsv	ALT_QUAL_test.1.tsv	ALT_FREQ_test.1.tsv	TOTAL_DP_test.1.tsv	PVAL_test.1.tsv	PASS_test.1.tsv	REF_DP_test.2.tsv	REF_RV_test.2.tsv	REF_QUAL_test.2.tsv	ALT_DP_test.2.tsv	ALT_RV_test.2.tsv	ALT_QUAL_test.2.tsv	ALT_FREQ_test.2.tsv	TOTAL_DP_test.2.tsv	PVAL_test.2.tsv	PASS_test.2.tsv	REF_DP_test.3.tsv	REF_RV_test.3.tsv	REF_QUAL_test.3.tsv	ALT_DP_test.3.tsv	ALT_RV_test.3.tsv	ALT_QUAL_test.3.tsv	ALT_FREQ_test.3.tsv	TOTAL_DP_test.3.tsv	PVAL_test.3.tsv	PASS_test.3.tsv	
test	139	T	A	id-test3	GCT	A	GCA	A	1	0	32	1	0	55	0.5	2	0.666667	FALSE	1	0	32	1	0	55	0.5	2	0.666667	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
test	320	A	T	NA	NA	NA	NA	NA	1	1	35	1	1	46	0.5	2	0.666667	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1	1	35	1	1	46	0.5	2	0.666667	FALSE
test	365	A	T	NA	NA	NA	NA	NA	0	0	0	1	1	27	1	1	1	FALSE	0	0	0	1	1	27	1	1	1	FALSE	0	0	0	1	1	27	1	1	1	FALSE
test	42	G	T	id-test4	CAG	Q	CAT	H	0	0	0	1	0	49	1	1	1	FALSE	0	0	0	1	0	49	1	1	1	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
test	42	G	T	id-testedit1	AGG	R	ATG	M	0	0	0	1	0	49	1	1	1	FALSE	0	0	0	1	0	49	1	1	1	FALSE	0	0	0	1	0	49	1	1	1	FALSE
test	69	T	G	id-testedit2	TTG	L	TGG	W	1	0	57	1	0	53	0.5	2	0.666667	FALSE	1	0	57	1	0	53	0.5	2	0.666667	FALSE	1	0	57	1	0	53	0.5	2	0.666667	FALSE
```

Description of fields

| No | Field                                             | Description                                              |
|:---|:--------------------------------------------------|----------------------------------------------------------|
|  1 | REGION                                            | Common region across all replicate variant tsv files     |
|  2 | POS                                               | Common position across all variant tsv files             |
|  3 | REF                                               | Common reference base across all variant tsv files       |
|  4 | ALT                                               | Common alternate base across all variant tsv files       |
|  5 | GFF\_FEATURE                                               | GFF feature used for the translation       |
|  6 | REF\_CODON                                               | The codon using the reference base       |
|  7 | REF\_AA                                               | Reference codon translated into amino acid       |
|  8 | ALT\_CODON                                               | Codon using the alternate base       |
|  9 | ALT\_AA                                               | Alternate codon translated into amino acid       |
|  10 | REF_DP_<rep1-tsv-file-name>                       | Depth of reference base in replicate 1                   |
|  11 | REF_RV_<rep1-tsv-file-name>                       | Depth of reference base on reverse reads in replicate 1  |
|  12 | REF_QUAL_<rep1-tsv-file-name>                     | Mean quality of reference base in replicate 1            |
|  13 | ALT_DP_<rep1-tsv-file-name>                       | Depth of alternate base in replicate 1                   |
|  14 | ALT_RV_<rep1-tsv-file-name>                       | Deapth of alternate base on reverse reads in replicate 1 |
|  15 | ALT_QUAL_<rep1-tsv-file-name>                     | Mean quality of alternate base in replicate 1            |
|  16 | ALT_FREQ_<rep1-tsv-file-name>                     | Frequency of alternate base in replicate 1               |
|  17 | TOTAL_DP_<rep1-tsv-file-name>                     | Total depth at position in replicate 1                   |
|  18 | PVAL_<rep1-tsv-file-name>                         | p-value of fisher's exact test in replicate 1            |
|  19 | PASS_<rep1-tsv-file-name>                         | Result of p-value <= 0.05 in replicate 1                  |
|  20 | Continue rows 10 - 19 for every replicate provided |                                                          |

Generate a consensus sequences from an aligned BAM file
----

To generate a consensus sequence iVar uses the output of `samtools mpileup` command. The mpileup output must be piped into `ivar consensus`. There are five parameters that can be set -  minimum quality(Default: 20), minimum frequency threshold(Default: 0), minimum depth to call a consensus(Default: 10), a flag to exclude nucleotides from regions with depth less than the minimum depth and a character to call in regions with coverage lower than the speicifed minimum depth(Default: 'N'). Minimum quality is the minimum quality of a base to be considered in calculations of variant frequencies at a given position. Minimum frequency threshold is the minimum frequency that a base must match to be called as the consensus base at a position. If one base is not enough to match a given frequency, then an ambigious nucleotide is called at that position. Minimum depth is the minimum required depth to call a consensus. If '-k' flag is set then these regions are not included in the consensus sequence. If '-k' is not set then by default, a 'N' is called in these regions. You can also specfy which character you want to add to the consensus to cover regions with depth less than the minimum depth. This can be done using -n option. It takes one of two values: '-' or 'N'.

As an example, consider a position with 6As, 3Ts and 1C. The table below shows the consensus nucleotide called at different frequencies.

| Minimum frequency threshold | Consensus |
|:----------------------------|-----------|
| 0 | A |
| 0.5 | A |
| 0.6 | A |
| 0.7 | W(A or T) |
| 0.9 | W (A or T) |
| 1 | H (A or T or C) |

If there are two nucleotides at the same frequency, both nucleotides are used to call an ambigious base as the consensus. As an example, consider a position wiht 6 Ts, 2As and 2 Gs. The table below shows the consensus nucleotide called at different frequencies.

| Minimum frequency threshold | Consensus |
|:----------------------------|-----------|
| 0 | T |
| 0.5 | T |
| 0.6 | T |
| 0.7 | D(A or T or G) |
| 0.9 | D(A or T or G) |
| 1 | D(A or T or G) |

The output of the command is a fasta file with the consensus sequence and a .txt file with the average quality of every base used to generate the consensus at each position. *For insertions, the quality is set to be the minimum quality threshold since mpileup doesn't give the quality of bases in insertions.*

Command:

```

ivar consensus

Usage: samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 

Note : samtools mpileup output must be piped into ivar consensus

Input Options    Description
           -q    Minimum quality score threshold to count base (Default: 20)
           -t    Minimum frequency threshold(0 - 1) to call consensus. (Default: 0)
                 Frequently used thresholds | Description
                 ---------------------------|------------
                                          0 | Majority or most common base
                                        0.2 | Bases that make up atleast 20% of the depth at a position
                                        0.5 | Strict or bases that make up atleast 50% of the depth at a position
                                        0.9 | Strict or bases that make up atleast 90% of the depth at a position
                                          1 | Identical or bases that make up 100% of the depth at a position. Will have highest ambiguities
           -m    Minimum depth to call consensus(Default: 10)
           -k    If '-k' flag is added, regions with depth less than minimum depth will not be added to the consensus sequence. Using '-k' will override any option specified using -n 
           -n    (N/-) Character to print in regions with less than minimum coverage(Default: N)

Output Options   Description
           -p    (Required) Prefix for the output fasta file and quality file
```

Example Usage:
```
samtools mpileup -d 1000 -A -Q 0 test.bam | ivar consensus -p test -q 20 -t 0
```

The command above will produce a test.fa fasta file with the consensus sequence and a test.qual.txt with the average quality of each base in the consensus sequence.

Get primers with mismatches to the reference sequence
----

iVar uses a .tsv file with variants to get the zero based indices(based on the BED file) of mismatched primers. This command requires another .tsv file with each line containing the left and right primer names separated by a tab. This is used to get both the primers for an amplicon with a single mismatched primer. The output is a text file with the zero based primer indices delimited by a space. The output is written to a a text file using the prefix provided.

Command:
```
ivar getmasked
Usage: ivar getmasked -i <input-filtered.tsv> -b <primers.bed> -f <primer_pairs.tsv> -p <prefix>
Note: This step is used only for amplicon-based sequencing.

Input Options    Description
           -i    (Required) Input filtered variants tsv generated from 'ivar filtervariants'
           -b    (Required) BED file with primer sequences and positions
           -f    (Required) Primer pair information file containing left and right primer names for the same amplicon separated by a tab
Output Options   Description
           -p    (Required) Prefix for the output text file

```

Example BED file

```
Puerto	28	52	400_1_out_L	60	+
Puerto	482	504	400_1_out_R	60	-
Puerto	359	381	400_2_out_L	60	+
Puerto	796	818	400_2_out_R	60	-
Puerto	658	680	400_3_out_L*	60	+
Puerto	1054	1076	400_3_out_R*	60	-
.
.
.
.
```

Example primer pair information file
```
400_1_out_L    400_1_out_R
400_2_out_L    400_2_out_R
400_3_out_L    400_3_out_R
.
.
.
.
```

Example Usage:
```
ivar getmasked -i test.filtered.tsv -b primers.bed -f pair_information.tsv -p test.masked.txt
```

The command above produces an output file - test.masked.txt.

Example Output:

```
1 2 7 8
```

Remove reads associated with mismatched primer indices
----

This command accepts an aligned and sorted BAM file trimmed using `ivar trim` and removes the reads corresponding to the supplied primer indices, which is the output of `ivar getmasked` command. Under the hood, `ivar trim` adds the zero based primer index(based on the BED file) to the BAM auxillary data for every read. Hence, ivar removereads will only work on BAM files that have been trimmed using `ivar trim`.

Command:
```
ivar removereads

Usage: ivar removereads -i <input.trimmed.bam> -p <prefix> -t <text-file-with-primer-indices>
Note: This step is used only for amplicon-based sequencing.

Input Options    Description
           -i    (Required) Input BAM file  trimmed with ivar trim. Must be sorted and indexed, which can be done using sort_index_bam.sh
           -t    (Required) Text file with primer indices separated by spaces. This is the output of getmasked command.

Output Options   Description
           -p    (Required) Prefix for the output filtered BAM file

```

Example Usage:
```
ivar trim -i test.bam -p test.trimmed
ivar removereads -i test.trimmed.bam -p test.trimmed.masked.bam -t test.masked.txt
```

The `ivar trim` command above trims test.bam and produced test.trimmed.bam with the primer indice data added. The `ivar removereads` command produces an output file - test.trimmed.masked.bam after removing all the reads corresponding to primer indices - 1,2,7 and 8.

(Experimental) trimadapter
----

**Note: This feature is under active development and not completely validated yet.**

trimadapter in iVar can be used to trim adapter sequences from fastq files using a supplied fasta file.
