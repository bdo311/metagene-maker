#!/usr/bin/python
# metagene_maker.py
# 8/29/14, last updated 12/23/14
# makes metagenes for bedgraphs and regions according to a configuration file

# SNF TODO 
# - to this point, report a histogram of region sizes for processed regions in region space (not chr space) 

import os, glob, csv, re, collections, math, multiprocessing, sys, random, subprocess, logging, argparse
import numpy as np
from datetime import datetime
from binning_functions import *
from merge_bins import *
csv.register_dialect("textdialect", delimiter='\t')

# parser
parser = argparse.ArgumentParser(description="metagene-maker: obtain average profiles of NGS datasets over your favorite regions", epilog="Example: python metagene_maker.py config.txt ChIP_exp outputDir/")
parser.add_argument('config_file', metavar='config_file', help='required configuration file')
parser.add_argument('prefix', metavar='prefix', help="Prefix of output files")
parser.add_argument('outputDir', metavar='output_directory', help="Directory where output folders will be written")
parser.add_argument('-l', metavar='binLength', type=int, help="Bases per window when processing bedgraph. Default is 2,000,000.", default=2000000)
parser.add_argument('-p', metavar='processors', type=int, help="Number of cores to use. Default is 4.", default=4)
args = parser.parse_args()
config_file = args.config_file
numProcs = args.p
name = args.prefix
binLength = args.l
parentDir = args.outputDir
if parentDir[0] != '/': parentDir = os.getcwd() + '/' + parentDir

# log file
logger=logging.getLogger('')
logger.setLevel(logging.INFO)
fh=logging.FileHandler('metagene.log')
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
fh.setFormatter(formatter)
ch=logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.addHandler(fh)

def readConfigFile(fn):
	# Input: a configuration file with folders and regions
	# Output: these parameters in a hashmap format

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	# folders
	folders = {}
	for row in reader:
		if len(row)==0: continue
		if 'regionType' in row[0]: break
		if '#' in row[0] or row[0]=='': continue
		folders[row[0]] = row[1:]

	# regions
	regions = {}
	for row in reader:
		if len(row)==0: continue
		if '#' in row[0] or row[0]=='': continue
		regions[row[0]] = row[1:]

	return folders, regions

def processFolders(parentDir, folders, regions):
	folderToGraph = {}
	allChroms = set()
	for folder in folders:
		# setting up folders
		os.chdir(parentDir)
		if not glob.glob(folder + '/'): os.mkdir(folder)
		os.chdir(folder)
		
		if not glob.glob('bins/'): os.system("mkdir bins")
		os.chdir('bins')
		#os.system("rmdir *") # removing empty directories
		for r in regions:
			if not glob.glob(r): os.system("mkdir " + r)

		# splitting up bedgraph if not done already
		os.chdir("..")
		if not glob.glob("bedGraphByChr/"): os.system("mkdir bedGraphByChr")
		
		os.chdir("bedGraphByChr")
		#os.system("rm -f *.bedGraph")
		logger.info("\nSplitting up bedgraph for %s", folder)
		cmd = "gawk '{print >> $1\".bedGraph\"}' " + folders[folder][0]
		logger.info(cmd)
		os.system(cmd)

		# adding chromosomes to list of all chroms
		files = [os.path.basename(fn) for fn in glob.glob('*.bedGraph')]
		chrs = [x.replace('.bedGraph', '') for x in files]
		allChroms.update(set(chrs))		
				
		# making folder to bedgraph relationship
		binFolder = parentDir + '/' + folder + '/bins/'
		graphFolder = parentDir + '/' + folder + '/bedGraphByChr/'
		folderToGraph[folder] = [binFolder, graphFolder, folders[folder][1]] #bin folder --> [graph folder, strand]

	return folderToGraph, allChroms

# checks whether file is a bed file
def isBed(row):
	if 'chr' not in row[0]: return False
	try:
		a,b=int(row[1]), int(row[2]) #start/end ok?
		
		if len(row)>11: 
			if ',' not in row[10] or ',' not in row[11]: return False #BED12s must have commas
			c,d =len(row[10].split(',')), len(row[11].split(',')) #number of blocks must be the same
			if c!=d: return False
	except: return False
	return True
	
def getChrToRegion(fn, header):
	with open(fn, 'r') as ifile:	
		regions = collections.defaultdict(lambda: []) # by chromosome
		if header: ifile.next()
		
		counter = 0
		for line in ifile:
			counter += 1
			row = line.split()
			if not isBed(row):
				logger.info("Line %d of file %s is not in BED format. Exiting.", counter, fn)
				exit()
			regions[row[0]].append(row)
			
	# is region BED6 or BED12?
	rowLen = len(regions[regions.keys()[0]][0])
	bedType = 'BED6' if rowLen<11 else 'BED12'
	
	return regions, bedType

def processRegions(regions):
	regionToChrMap = {}
	regionToBedType = {}
	for region in regions:
		info = regions[region]
		loc = info[0]
		isHeader = True if info[1] == 'y' else False
		
		# check that number of bins and the extension size will be handled gracefully
		numBins = int(info[4])
		extendRegion = info[5]
		sideExtension = int(info[6])

		if sideExtension and extendRegion=='y': 
			logger.info("Region %s cannot be extended two different ways. Exiting.", region)
			exit()
		if sideExtension:
			if numBins % 4 != 0: 
				logger.info("Number of bins in region %s must be a multiple of 4 when using fixed side extensions (column sideExtensions). Exiting.", region)
				exit()
			if sideExtension % (numBins/4) != 0: 
				logger.info("Size of fixed extension in region %s (column sideExtensions) must be a multiple of the number of bins divided by 4. Exiting.", region)
				exit()
		regionInfo = getChrToRegion(loc, isHeader)
		regionToChrMap[region] = regionInfo[0]
		regionToBedType[region] = regionInfo[1]
		if regionInfo[1] == 'BED12':
			if extendRegion == 'y': 
				logger.info("Spliced regions in %s cannot be extended using extendRegion. Exiting.", region)
				exit()
	return regionToChrMap, regionToBedType


def main():			
	# read config file
	folders, regions = readConfigFile(config_file)
	logger.info("\nRead configuration file")
	
	# processing folders and bedgraphs
	if not glob.glob(parentDir): os.system("mkdir " + parentDir)
	folderToGraph, allChroms = processFolders(parentDir, folders, regions)
	allChroms = list(allChroms)
	logger.info("\nProcessed folders: %s", ', '.join(folderToGraph.keys()))
	for f in folderToGraph:
		logger.info('%s: %s', f, folderToGraph[f][0])
	
	# processing regions, checking that they are valid bed files
	regionToChrMap, regionToBedType = processRegions(regions)
	logger.info("\nProcessed regions: %s", ', '.join(regionToChrMap.keys()))
	for r in regionToChrMap:
		logger.info('%s: %s %s', r, regionToBedType[r], regions[r][0])

	# making bins
	logger.info("\nReading in bedgraphs and making profiles for each region")
	
	for folder in folderToGraph:
		# if my bedgraph is stranded and my regions are stranded, only 
		# use the regions that correspond to the bedgraph strand
		[binFolder, graphFolder, folderStrand] = folderToGraph[folder]
		for i in range(len(allChroms)):
			chroms = allChroms[(numProcs*i):(numProcs*(i+1))]
			if chroms == []: break
			procs=[]
			for chrom in chroms:
				p = multiprocessing.Process(target=processEachChrom, args=(chrom, binFolder, graphFolder, folderStrand, binLength, regions, regionToChrMap, regionToBedType))
				p.start()
				procs.append(p)		
			for proc in procs: proc.join()		

	exit()
	# merging bins for each chromosome, then make metagene
	logger.info("\nMaking metagenes")
	folders = folderToGraph.keys() 
	numPerProc = len(folders)/numProcs + 1 if len(folders) > numProcs else 1
	numJobs = numProcs if (numProcs < len(folders)) else len(folders)
	procs = []

	for i in range(numJobs): 
		p = multiprocessing.Process(target=folderWorker, args=(i * numPerProc, (i + 1) * numPerProc, folders, folderToGraph, regions))
		procs.append(p)
		p.start()
	for p in procs: p.join()

	# merging all files, and writing average files
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	os.chdir(parentDir)
	if not glob.glob("averages"): os.system("mkdir averages")
	for region in regions:
		for folder in folderToGraph:
			binFolder = folderToGraph[folder][0]
			os.chdir(binFolder + '/' + region + '/')
			fn = "avgraw_" + folder + "_" + region 
			regionToFolderAvgs[region][folder] = processFile(fn)
		writeFile(name + '_' + region, regionToFolderAvgs[region], parentDir + '/averages/')

	logger.info("\nDone!\n")
	
if __name__ == '__main__':
	main()

