# CLOMP - CLinically Okay Metagenomic Pipeline
![CLOMP](https://github.com/FredHutch/CLOMP/workflows/CLOMP/badge.svg?branch=local_ngs_server)

A entirely open source, fast, accurate, multi-sample and highly configurable end to end metagenomics pipeline for the people.

# Introduction

This is the publicly available source code and documentation for CLOMP - UW Virology's fully functional metagenomic pipeline. CLOMP takes raw sequence read files and taxonomically assigns as many reads as possible. We have attempted to streamline the setup for both local and cloud use as much as possible by writing the pipeline in `Nextflow`.

Broadly the execution of the pipeline is broken down into four steps:
### 1. Read quality filtering and adapter trimming 

Adapter adn quality trimming is done using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) with a pre-made comprehensive Illumina adapter file. The parameters for trimming are adjustable via the  `--TRIMMOMATIC_OPTIONS` flag. 

### 2. Host read subtraction

Trimmed reads are aligned to a host genome using [bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3). Aligned reads are removed from further analysis in the interest of speed and computational cost. This pipeline is primarily designed for human clinical samples so by default we align to the Human Genome (HG38). Depending on the sample type, this may remove almost all reads. To use a different host genome, provide a path to a folder containing the built `Bowtie2` index files using the flag `--BWT_DB_LOCATION`, and the `Bowtie2` prefix for the files with `--BWT_DB_PREFIX `

### 3. Alignment to NT

All remaining reads are aligned to the NCBI NT database using the [SNAP](http://snap.cs.berkeley.edu/) aligner, chosen for its speed and accuracy.`SNAP` alignment options can be modified using the `--SNAP_OPTIONS` flag. Details on building the SNAP database are available below. 

### 4. Tiebreaking and visualization

This step taxonomically classifies each read based on the previous `SNAP` alignment to the NT database. The custom tiebreaking logic amalgamates all taxa the read aligned to, and assigns it to the lowesst common ancestor possible. This may be at any level of phylogeny (species, genus, etc). The output `.tsv` files can be visualized using [Pavian](https://github.com/fbreitwieser/pavian).

# Installation
## Required System Specs
This pipeline can run on any operating system capable of running `Nextflow`. However, in our experience the hardware specifications required are largely determined by the size of your `SNAP` databases. To use our pre-built databases, at a minimum 2.5tb of disk space and 244gb of ram are required. 

## Installation and setup
`CLOMP` was built with portability and modularity in mind. As such, all that is required is the latest version of `Nextflow` to run the pipeline, and `Docker` to fetch the required images. If you would like to use our databases, `AWScli` is also required to download those. Please note that If using our databases, the first run will be extremely slow as the databases are large and will require time to download. 

2. Each SNAP database has a different composition we go through all assignments and throw out any that are more than 6 edit distance away from the best hit. This essentially just throws out the case where we get one pretty tenuous hit in one database, this is much more common for less heavily sequenced organisims. 

3. We then take the most specific taxonomical assignment that has 100% agreement. 

The other tiebreaking options that are currently included are `90` and `oneoff`. 
 * 90: The most specific taxid shared between ~90% of the remaining alignments is returned
 * oneoff: If there is not more than one other taxid assigned, call read most common taxid, otherwise intersect and pick most specific shared taxid.  Essentially, one weird/bad assignment cannot define/determine the tree intersection.

### Test data run script 
Test data is provided in the `test_data` folder to confirm the pipeline is installed and runnable. The following script should run when executed from the root of the `CLOMP` repository folder: 

     NXF_ANSI_LOG=false NXF_VER=20.01.0 nextflow run main.nf \
     -profile testing \
     --INPUT_FOLDER test_data/ \
     --INPUT_SUFFIX .fastq.gz  \
     --ENTREZ_EMAIL vpeddu@uw.edu \
     --OUTDIR test_data/out/  \
     -work-dir test_data/work/ \
     --SNAP_INDEXES "test_data/test_db_1/,test_data/test_db_2/,test_data/test_db_3/" \
     --TRIMMOMATIC_OPTIONS ":2:30:10 SLIDINGWINDOW:4:20"
     --TRIMMOMATIC_JAR_PATH test_data/trimmomatic-0.38.jar \
     -with-docker ubuntu:18.04 \
     --KRAKEN_DB_PATH test_data/kraken_db/ 
     --BWT_DB_LOCATION test_data/bt2/ \
     --TRIMMOMATIC_ADAPTER_PATH test_data/adapters.fa \
     --BWT_DB_PREFIX human_ref


### All arguments: 

#### Required Arguments:
      * --INPUT_FOLDER Folder containing input files (FASTQ format)
      * --ENTREZ_EMAIL Email adress requried for interfacing with NCBI servers. Please provide your own email address.
      * --OUTDIR The output directory where the results will be saved
#### Optional Arguments:
      * --PAIRED_END If specified, look for paired-end files with { _R1_ / _R2_ }
      * --INPUT_SUFFIX Suffix of input files (default: .fastq.gz)
      * --SNAP_INDEXES Comma-delimited list of SNAP indexed references
      * --SNAP_OPTIONS Options used to run SNAP (must enclose in quotes). Default: -mrl 65 -d 9 -h 30000 -om 1 -omax 20 -map
      * --SNAP_BATCHSIZE Number of samples to align in parallel over each SNAP index shard, in a batch
      * --HOST_FILTER_FIRST If specified, perform host filtering prior to trimming
      * --SECOND_PASS Should we do a more sensitive second pass of host filtering? (default: false)
      * --TRIMMOMATIC_OPTIONS Options used to run Trimmomatic (default: ':2:30:10 HEADCROP:10 SLIDINGWINDOW:4:20 CROP:65 MINLEN:65')
      * --TRIMMOMATIC_JAR_PATH Path to Trimmomatic executable (default: "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/trimmomatic-0.38.jar")
      * --TRIMMOMATIC_ADAPTER_PATH Path to Trimmomatic adapters (default: "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/adapters.fa")
      * --BWT_DB_PREFIX Prefix for the files in the host genome directory (default: hg38)
      * --BWT_DB_LOCATION Host genome directory (default: s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/hg38/)
      * --BWT_SECOND_PASS_OPTIONS Bowtie2 options used in the second pass analysis (default: '-D 35 -R 5 -N 1 -L 19 -i S,1,0.50')
      * --KRAKEN_DB_PATH The directory containing the Kraken database (default: s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/kraken_db/)
      * --LOGIC This controls the different types of tie-breaking logic that are baked in to CLOMP: strict/90/oneoff (default: strict)
      * --BLAST_EVAL E-value cuttoff for BLAST (default: 0.001)
      * --BLAST_CHECK Check results with BLAST (default: false)
      * --BLAST_CHECK_DB Database used to check results with BLAST (default: None)
      * --FILTER_LIST List of taxids to filter from tiebreaking (default: "[12908,28384,48479]")
      * --KRAKEN_DB_PATH Path to Kraken database (default: "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/kraken_db/")
      * --H_STRICT Assign reads with any host alignment to host (default: false)
      * --H_TAXID Host taxid (default: 9606)
      * --INCLUSION_TAXID Inclusion TAXID will pull all reads that get classified below this level but not downstream of the exclusion taxid (default: 2759)
      * --EXCLUSION_TAXID See above (default: 9604)
      * --BASE_DELIMITER Each sample needs a unique name and everything BEFORE the FIRST instance of this string will be the sample base name (default: '_')
      * --MIN_READ_CUTOFF The minimum number of reads assigned to a species level taxid that are required to build the sam file (default: 10)
      * --SAM_NO_BUILD_LIST Taxids that contain these taxids in their lineage will not be built (default: "[2759,77133]")
      * --ADD_HOST_FILTERED_TO_REPORT Should we add the number of reads called to the host to the final Pavian output? (default: true)
      * --HOST_FILTER_TAXID The taxid for which we should classify host filtered reads as when adding them to the Pavian report (default: 9606)
      * --WRITE_UNIQUES Should program write a new folder named after the sample and write all unique reads that mapped to each taxid into this folder? VERY useful for debugging tie-breaking (default: true)
      * --EDIT_DISTANCE_OFFSET The edit distance offset is the maximum difference in edit distance we will accept for alignments (default: 6)
      * --BUILD_SAMS Should we build SAM files aligning reads to species level assignments? (default: false)
### Ask me about the license ###

I believe in free software. I think it would be extrememly disingenuous to the bioinformatics community and the greater scientific public to allow what is essentially a parsing script for multiple sam files to be put under a restrictive license. This entire software is based off publically available databases and tools - I also think that this pipeline (as in the order of the steps and the arguments) works quite well and could be widely used, hopefully making people's lives easier. 
