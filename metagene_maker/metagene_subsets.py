# metagene_subsets.py
# takes already run metagene maker and applies user generated subsets
# 1/30/15

import os, glob, csv, re, collections, math, multiprocessing, sys, random, subprocess, logging, argparse
from datetime import datetime
csv.register_dialect("textdialect", delimiter='\t')

# parser
parser = argparse.ArgumentParser(description="metagene-subsets: obtain average profiles of subsets of regions already processed by metagene maker", epilog="Example: python metagene_subsets.py config.txt ChIP_exp outputDir/")
parser.add_argument('config_file', metavar='config_file', help='required configuration file')
parser.add_argument('prefix', metavar='prefix', help="Prefix of output files")
parser.add_argument('outputDir', metavar='output_directory', help="Directory where output folders will be written")
parser.add_argument('-p', metavar='processors', type=int, help="Number of cores to use. Default is 4.", default=4)
args = parser.parse_args()
config_file = args.config_file
numProcs = args.p
prefix = args.prefix
parentDir = args.outputDir
if parentDir[0] != '/': parentDir = os.getcwd() + '/' + parentDir


def readConfigFile(fn):
	# Input: a configuration file with folders and regions
	# Output: these parameters in a hashmap format

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	# folders
	# stranded bedgraphs must come in pairs. all pairs must have a plus and a minus
	folders = {}
	for row in reader:
		if len(row)==0: continue
		if 'regionType' in row[0]: break
		if '#' in row[0] or row[0]=='': continue
		
		strand = row[2]
		assert(strand == '0' or strand == '+' or strand == '-')

		folderName = row[0]
		folderPath = parentDir + folderName
		if not glob.glob(folderPath):
			logger.info("Folder %s not found in parent directory %s. Exiting.", folderName, parentDir)
			exit()
		folders[row[0]] = [folderPath, strand]		

	# regions
	regions = {}
	for row in reader:
		if len(row)==0: continue
		if '#' in row[0] or row[0]=='': continue
		regions[row[0]] = row[1:]

	return folders, regions
	
def getColumnMean(dir, isMinus):
	a=pd.read_table(dir + "allchr_sorted.txt", header=None)
	b=a[range(7,len(a.columns) + 1]
	x=b.mean(axis=0)
	if isMinus: return list(-x)
	return list(x)
	
def main():
	# 1. read config file
	folders, regions = readConfigFile(config_file)
	
	# 2. make subset folders and allchr.txt
	for region in regions:
		logger.info("Processing %s", region)
		subsetList = region[0]
		origRegion = region[1]
		
		# get list of region names to subset
		names = set()
		with open(subsetList, 'r') as ifile:
			for line in ifile:
				names.add(line.rstrip())
		
		for folder in folders:
			logger.info("\t%s", folder)
			folderPath = folders[folder]
			origRegionFolder = folderPath + '/bins/' + origRegion
			if not glob.glob(origRegionFolder):
				logger.info("Region type %s not found in sample %s. Exiting.", folder, region)
				exit()
			
			allChrSorted = origRegionFolder + "/allchr_sorted.txt"
			assert(glob.glob(allChrSorted))
			
			# make new bin folder
			outputFolder = folderPath + '/bins/' + region
			os.mkdir(outputFolder)
			os.chdir(outputFolder)
			
			# subset
			with open("allchr.txt", 'w') as ofile, open(allChrSorted, 'r') as ifile:
				reader = csv.writer(ifile, 'textdialect')
				writer = csv.writer(ofile, 'textdialect')
				
				for row in reader:
					rowName = row[3]
					if rowName in names: writer.writerow(row)
					
			# sort and remove
			os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")
			os.system("rm -f allchr.txt")
	
	# 3. getting metagenes and writing average files.
	logger.info("\nWriting average files")
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	os.chdir(parentDir)
	for region in regions:
		for folder in folderToGraph:
			binFolder = folderToGraph[folder][0] + '/bins/'
			isMinus = (folderToGraph[folder][1] == '-')
			dir = binFolder + '/' + region + '/'
			regionToFolderAvgs[region][folder] = getColumnMean(dir, isMinus)
		logger.info("%s_%s", prefix, region)
		writeFile(prefix + '_' + region, regionToFolderAvgs[region], parentDir + '/averages/')

	logger.info("\nDone!\n")
	

if __name__ == '__main__':
	main()


