metagene-maker
==============

Makes metagene plots for any bedgraph over given regions. Python and R should be installed.

Simple start
----------

1. Clone project to new directory
2. Make config file (see below)
3. Run: `metagene_maker.py <config file>`, either in `screen` or `nohup`
4. Output: excel spreadsheets for each region in a new `averages` folder in the user-provided parent directory

Configuration file
--------

### Parameters
**name:** the name appended to all output files

**parentDir:** the folder where all input and output subfolders are located

**organism:** either hg19 or mm9 depending on organism

**threads:** number of processors used

### Bedgraph columns
**folder:** the name of the sample (should also be the name of the folder where sample-specific intermediate files will be made)

**bedGraphLoc:** path to bedgraph

**stranded:** + if plus only, - if minus only, 0 if no strand information

### Region columns
**regionType:** name of region

**fileLoc:** path to file specifying the regions of interest

**header:** y if header, n if no header

**chrCol:** 0-indexed column number of chromosome designations

**nameCol:** 0-indexed column number of name designations

**startCol:** 0-indexed column number of start designations

**stopCol:** 0-indexed column number of stop designations

**stranded:** y if directional (i.e. TSS's), n if not directional (i.e. enhancers)

**strandCol:** 0-indexed column number of strand designations; put 0 if no strand

**numCols:** number of columns in file

**limitSize:** y if only regions >200bp and <200kb should be considered; n if no limitation

**numBins:** number of bins. anywhere between 100 and 1000 is good

**extendRegion:** y if regions in the bed file should be extended 1x upstream and downstream; n otherwise


