metagene-maker
==============

Makes metagene plots for bedgraphs over given regions in bed files for any organism. Regions can be continuous or spliced. Useful in analysis of ChIRP-seq, ChIP-seq, GRO-seq, ATAC-seq, iCLIP, ribosome profiling, RNA-seq, and other NGS datasets.

**Table of Contents** 

- [metagene-maker](#)
	- [Installation](#)
	- [Usage](#)
	- [Dependencies](#)
	- [Making BED files](#)
	- [Configuration file (see example.conf for an example)](#)
		- [Bedgraph columns](#)
		- [Region columns](#)
	- [Output directory structure](#)
		- [Parent directory](#)
		- [Subfolders](#)
		- [Contents of each sample folder](#)

Installation
----------

1. Go to 'releases' above and download the latest tar.gz file. Unzip with `tar xvzf metagene-maker-0.x.tar.gz>`
2. Alternatively, you can clone this git repository using `git clone`.
3. Go into the folder: `cd <metagene-maker-0.x>`
4. Make sure you have the needed dependencies (below). Install: `sudo python setup.py install`. If you do not have sudo privileges, run `python setup.py install --user` or `python setup.py install --prefix=<desired directory>`. Be sure that the python you use to run `setup.py` is version 2.7; scripts WILL NOT WORK with lower versions (2.4, 2.5).

Usage
-----
1. Make config file (see below)
2. Ensure that you have a bedgraph for every sample you want to analyze.
3. Ensure that you have properly formatted BED6/12 files for every region for which you want to build average profiles. You can make these with the included `extractTranscriptRegions` module (see below).
4. Run: `metagene_maker <config file> <name> <outputDir>` where <config file> is the configuration file you make using `example.conf` (provided) as the template. Instructions for making configuration file are below. Run this either in `screen` or `nohup`.
5. Output: tab delimited files for each region in a new `averages` folder in the user-provided output directory, as well as raw files named `allchr_sorted.txt` in each subfolder that contains binned profiles for each region and can be used for custom analysis.


usage: `metagene_maker [-h] [-l binLength] [-p processors] config_file prefix output_directory`

example: `metagene_maker -p 10 -l 500000 config/test.txt M3_ChIP chip/`

positional arguments: | explanation
--------------------|----------------------------
  config_file    |   required configuration file
  prefix         |   Prefix of output files
  output_directory |  Directory where output folders will be written

optional arguments: | explanation
-------------------|-------------------------------------------
  -h, --help      |  show this help message and exit
  -l binLength    |  Bases per window when processing bedgraph. Default is 2,000,000.
  -p processors   |  Number of cores to use. Default is 4.


Dependencies
--------

1. Python (>=2.7)
2. Numpy (a python module) (>=1.7)
3. Pandas (a python module) (>=0.14)

At least 4 GB RAM if your largest bedgraph is 1 GB and you use 4 cores (empirical rule: n cores * m GB bedgraph --> mn GB RAM needed)

Making BED files
--------

You can supply your own BED6/12 files or use genome-wide BED files made using an included script, extractTranscriptRegions. You can start from either GTF files or files downloaded from UCSC as follows:

### Instructions for GTF files (use --gtf flag):
1. Download the GTF file into the desired directory.

### Instructions for UCSC files (use --ucsc flag):

1. From UCSC Genome Browser, go to Table Browser and choose your favorite organism/assembly. Choose "Genes and Gene Predictions" in 'group' and one of the gene tracks (we recommend UCSC Genes, Ensembl, or RefSeq).
2. Choose 'selected fields from primary and related tables' for 'output format'. 
3. Columns MUST be in this format: 
    - name
    - chrom
    - strand (+/-)
    - txStart
    - txEnd
    - cdsStart
    - cdsEnd
    - exonCount
    - exonStarts
    - exonEnds 
    - score
    - name2
4. Download the file.

### Running the script

Run `extractTranscriptRegions -i <gene_file.txt> -o <output_prefix> [--ucsc|--gtf]`. Output will be a list of bed files for UTRs, CDS's, exons, introns, splice sites, TSS's, and TES's that can be used for metagene-maker.

Configuration file (see example.conf for an example)
--------

### Bedgraph columns
**folder:** the name of the sample (should also be the name of the folder where sample-specific intermediate files will be made)

**bedGraphLoc:** absolute path to bedgraph

**stranded:** + if plus only, - if minus only, 0 if no strand information. IMPORTANT if your regions are also strand specific.

**pairName:** If a bedgraph is stranded, it must be part of a pair of bedgraphs (one + and one -) that share the same pairName. 

### Region columns

**regionType:** name of region

**fileLoc:** absolute path to file specifying the regions of interest

**limitSize:** y if only regions >200bp and <200kb should be considered; n if no limitation

**numBins:** number of bins. Use 1 to get the average coverage across the entire region. To make plots, anywhere between 100 and 500 is sufficient.

**extendRegion:** y if regions in the bed file should be extended 1x upstream and downstream; n otherwise

**sideExtension:** number of nt's to extend on each side of the provided regions. Default is 0.

**sideNumBins:** number of bins to allocate for the side extensions (this number must be less than half of numBins

Output directory structure
------

### Parent directory
the directory specified in the configuration file; contains all files generated by this pipeline

### Subfolders

**averages**: contains one file for each region type; each file is an Excel spreadsheet with graphable metagenes for each sample

**\<sample\>**: a folder for each sample, named as described in the config file; contains intermediate files described below

### Contents of each sample folder

**bedGraphByChr:** bedgraphs split by chromosome

**bins:** in this folder, there are subfolders for each region type, containing profiles for each instance of the region, intermediate RData files, and metagene plots

