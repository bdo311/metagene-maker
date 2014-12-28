#!/usr/bin/python2.7
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
prefix = args.prefix
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
	# stranded bedgraphs must come in pairs. all pairs must have a plus and a minus
	folders = {}
	folderPairs = collections.defaultdict(lambda: ['', ''])
	for row in reader:
		if len(row)==0: continue
		if 'regionType' in row[0]: break
		if '#' in row[0] or row[0]=='': continue
		folders[row[0]] = row[1:]
		
		strand = row[2]
		if strand == '0': continue #unstranded
		if len(row) <= 3 or row[3] == '':
			logger.info("Bedgraph %s must be part of a pair because it is stranded. Exiting.", row[0])
			exit()
		if strand == '+': folderPairs[row[3]][0] = row[0]
		elif strand == '-': folderPairs[row[3]][1] = row[0]
		else:
			logger.info("The strand for bedgraph %s must be '+' or '-' not %s. Exiting.", row[0], row[2])
			exit()
		
	# make sure each pair has a + and a - bedgraph
	for pairName in folderPairs:
		pair = folderPairs[pairName]
		if pair[0] == '' or pair[1] == '':
			logger.info("Pair %s must have both a '+' and a '-' bedgraph. Exiting.", pairName)
			exit()

	# regions
	regions = {}
	for row in reader:
		if len(row)==0: continue
		if '#' in row[0] or row[0]=='': continue
		regions[row[0]] = row[1:]

	return folders, folderPairs, regions

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
		logger.info("Splitting up bedgraph for %s", folder)
		if not glob.glob("done"): 
			os.system("rm -f *.bedGraph")
			cmd = "gawk '{print >> $1\".bedGraph\"}' " + folders[folder][0]
			logger.info(cmd)
			os.system(cmd)
			os.system("touch done")

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
		
		if row[5] != '+' and row[5] != '-': return False #must have strand info
		
		if len(row)>11: 
			if ',' not in row[10] or ',' not in row[11]: return False #BED12s must have commas
			c,d =len(row[10].split(',')), len(row[11].split(',')) #number of blocks must be the same
			if c!=d: return False
	except: return False
	return True
	
def getChrToRegion(fn):
	os.system("perl -p -i -e \"s/\r\n/\n/g\" " + fn) #replace newlines with the right newline
	with open(fn, 'r') as ifile:	
		regions = collections.defaultdict(lambda: []) # by chromosome
		
		counter = 0
		for line in ifile:
			counter += 1
			row = line.split()
			if not isBed(row):
				if counter == 1: continue #header
				logger.info("Line %d of file %s is not in BED6/12 format. Exiting.", counter, fn)
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
		
		# check that number of bins and the extension size will be handled gracefully
		numBins = int(info[2])
		extendRegion = info[3]
		sideExtension = int(info[4])
		sideNumBins = int(info[5])

		if sideExtension and extendRegion=='y': 
			logger.info("Region %s cannot be extended two different ways. Exiting.", region)
			exit()
		if sideExtension:
			if sideNumBins >= numBins/2: 
				logger.info("Number of bins for each side extension (%d) in region %s must be less than half the total number of bins (%d). Exiting.", sideNumBins, region, sideNumBins)

		regionInfo = getChrToRegion(loc)
		regionToChrMap[region] = regionInfo[0]
		regionToBedType[region] = regionInfo[1]
		if regionInfo[1] == 'BED12':
			if extendRegion == 'y': 
				logger.info("Spliced regions in %s cannot be extended using extendRegion. Exiting.", region)
				exit()
	return regionToChrMap, regionToBedType


def main():
	# 1. read config file
	folders, folderPairs, regions = readConfigFile(config_file)
	
	# 2. processing folders and bedgraphs
	logger.info("\nProcessing bedgraphs")
	if not glob.glob(parentDir): os.system("mkdir " + parentDir)
	folderToGraph, allChroms = processFolders(parentDir, folders, regions)
	allChroms = list(allChroms)
	logger.info("\nProcessed folders: %s", ', '.join(folderToGraph.keys()))
	for f in folderToGraph:
		logger.info('%s: %s', f, folderToGraph[f][0])
	
	# 3. processing regions, checking that they are valid bed files
	regionToChrMap, regionToBedType = processRegions(regions)
	logger.info("\nProcessed regions: %s", ', '.join(regionToChrMap.keys()))
	for r in regionToChrMap:
		logger.info('%s: %s %s', r, regionToBedType[r], regions[r][0])

	# 4. making bins
	logger.info("\nReading in bedgraphs and making profiles for each region")	
	for folder in folderToGraph:
		[binFolder, graphFolder, folderStrand] = folderToGraph[folder]
		for i in range(len(allChroms)):
			chroms = allChroms[(numProcs*i):(numProcs*(i+1))]
			if chroms == []: break
			procs=[]
			for chrom in chroms:
				p = multiprocessing.Process(target=processEachChrom, args=(chrom, binFolder, graphFolder, binLength, regions, regionToChrMap, regionToBedType))
				p.start()
				procs.append(p)		
			for proc in procs: proc.join()		

	# 5. making allchr_sorted.txt for each region for each dataset
	logger.info("\nMaking metagenes")	
	folders = folderToGraph.keys() 
	numPerProc = len(folders)/numProcs + 1 if len(folders) > numProcs else 1
	numJobs = numProcs if (numProcs < len(folders)) else len(folders)
	procs = []

	for i in range(numJobs): 
		p = multiprocessing.Process(target=concatChrs, args=(i * numPerProc, (i + 1) * numPerProc, folders, folderToGraph, regions))
		procs.append(p)
		p.start()
	for p in procs: p.join()
	
	# 6. combining stranded bedgraphs 
	logger.info("\nCombining stranded bedgraphs")
	procs = []
	for pair in folderPairs: #will be zero pairs if there are no stranded bedgraphs
		dir1 = parentDir + "/" + pair + "_sense/"
		dir2 = parentDir + "/" + pair + "_antisense/"
		if not glob.glob(dir1): os.system(" ".join(["mkdir", dir1, dir2]))
		os.chdir(dir1)
		if not glob.glob("bins"): os.system("mkdir bins")
		os.chdir(dir2)
		if not glob.glob("bins"): os.system("mkdir bins")
		
		# all regions will be broken up into sense and antisense. for non-stranded regions, treat as if it is (+)
		p = multiprocessing.Process(target = processPaired, args = (pair, folderPairs, regions, folderToGraph, parentDir))
		procs.append(p)
		p.start()
	for p in procs: p.join()
		
	# 7. remove original stranded bedgraph folders, and replace with sense/antisense	
	newFolderToGraph = {}
	for folder in folderToGraph: # get non stranded
		if folderToGraph[folder][2] == '0': newFolderToGraph[folder] = folderToGraph[folder]
	for pair in folderPairs: #make new entries for new folders
		senseName = pair + "_sense"
		binFolder = parentDir + '/' + senseName + "/bins/"
		graphFolder = parentDir + '/' + senseName + '/bins/'
		strand = '+'
		newFolderToGraph[senseName] = [binFolder, graphFolder, strand]
		
		asName = pair + "_antisense"
		binFolder = parentDir + '/' + asName + "/bins/"
		graphFolder = parentDir + '/' + asName + '/bins/'
		strand = '-'
		newFolderToGraph[asName] = [binFolder, graphFolder, strand]		
	folderToGraph = newFolderToGraph
	
	# 8. running R to get metagenes	
	folders = folderToGraph.keys() 
	numPerProc = len(folders)/numProcs + 1 if len(folders) > numProcs else 1
	numJobs = numProcs if (numProcs < len(folders)) else len(folders)
	procs = []

	for i in range(numJobs): 
		p = multiprocessing.Process(target=runRScript, args=(i * numPerProc, (i + 1) * numPerProc, folders, folderToGraph, regions))
		procs.append(p)
		p.start()
	for p in procs: p.join()

	# 9. merging all files, and writing average files. all antisense folders should have negative tracks
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	os.chdir(parentDir)
	if not glob.glob("averages"): os.system("mkdir averages")
	for region in regions:
		for folder in folderToGraph:
			binFolder = folderToGraph[folder][0]
			isMinus = (folderToGraph[folder][2] == '-')
			os.chdir(binFolder + '/' + region + '/')
			fn = "avgraw_" + folder + "_" + region 
			regionToFolderAvgs[region][folder] = processFile(fn, isMinus)
		writeFile(prefix + '_' + region, regionToFolderAvgs[region], parentDir + '/averages/')

	logger.info("\nDone!\n")
	
if __name__ == '__main__':
	main()

