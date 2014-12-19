metagene-maker
==============

Makes metagene plots for bedgraphs over given regions in bed files for human and mouse only. Useful in analysis of ChIRP-seq, ChIP-seq, GRO-seq, and other NGS datasets.

Simple start
----------

1. Clone project to new directory by executing `git clone https://github.com/bdo311/metagene-maker.git`. No need to install anything.
2. Make config file (see below)
3. Ensure that you have a bedgraph for every sample you want to analyze.
4. Ensure that you have properly formatted BED6 files for every region for which you want to build average profiles.
5. Run: `python metagene_maker.py <config file>` where <config file> is the configuration file you make using `example.conf` (provided) as the template. Instructions for making configuration file are below. Run this either in `screen` or `nohup`.
6. Output: tab delimited files for each region in a new `averages` folder in the user-provided parent directory

Dependencies
--------

Base Python (>2.7) and R (>3.0)

Rscript should be callable from the command line

At least 4 GB RAM if your largest bedgraph is 1 GB and you use 4 cores (empirical rule: n cores * m GB bedgraph --> mn GB RAM needed)

Todo
--------

Still doesn't work too well for stranded bedgraphs (i.e. GRO-seq data). Implementation is in progress.

Want to make metagenes that concatenate previously made metagenes (i.e. promoter, CDS, TES). Currently implemented as part of R script but not part of package yet.

Want to make metagenes for mRNAs (5'UTR, CDS, 3'UTR). Introns need to be thrown out.

Report a histogram of region sizes for processed regions in region space (not chr space)

Parse blocks for multi exon regions in the input bed file and turn these into a new object that has a method that can map bin space onto chr space and vice versa (SNF working on currently)


Configuration file (see example.conf for an example)
--------

### Parameters
**name:** the name appended to all output files

**parentDir:** the folder where all output subfolders will be located (see below for directory structure)

**organism:** either hg19 or mm9 depending on organism

**threads:** number of processors used

### Bedgraph columns
**folder:** the name of the sample (should also be the name of the folder where sample-specific intermediate files will be made)

**bedGraphLoc:** absolute path to bedgraph

**stranded:** + if plus only, - if minus only, 0 if no strand information. IMPORTANT if your regions are also strand specific.

### Region columns
**regionType:** name of region

**fileLoc:** absolute path to file specifying the regions of interest

**header:** y if header, n if no header

**stranded:** y if directional (i.e. TSS's), n if not directional (i.e. enhancers). IMPORTANT because some region profiles (like transcription start sites) have assymetrical shapes.

**limitSize:** y if only regions >200bp and <200kb should be considered; n if no limitation

**numBins:** number of bins. anywhere between 100 and 1000 is good

**extendRegion:** y if regions in the bed file should be extended 1x upstream and downstream; n otherwise

Directory structure
------

### Parent directory
the directory specified in the configuration file; contains all files generated by this pipeline

### Subfolders

**averages**: contains one file for each region type; each file is an Excel spreadsheet with graphable metagenes for each sample

**\<sample\>**: a folder for each sample, named as described in the config file; contains intermediate files described below

### Contents of each sample folder

**bedGraphByChr:** bedgraphs split by chromosome

**bins:** in this folder, there are subfolders for each region type, containing profiles for each instance of the region, intermediate RData files, and metagene plots

