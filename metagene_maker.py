#!/usr/bin/python
# metagene_maker.py
# 8/29/14
# makes metagenes for bedgraphs and regions according to a configuration file

# SNF TODO 
# - does not handle short (especially 1 nt long) regions - dies inside Rscript.  Use some cutoff value to trim super short things?  These may distort analysis 
# - to this point, report a histogram of region sizes for processed regions in region space (not chr space) 
# - parse blocks for multi exon regions in the input bed file and turn these into a new object that has a method that can map bin space onto chr space and vice versa 


import os, glob, csv, re, collections, math, multiprocessing, sys, random, subprocess
from datetime import datetime
from binning_functions import *
from merge_bins import *
csv.register_dialect("textdialect", delimiter='\t')

def readConfigFile(fn):
	# Input: a configuration file with parameters, folders, and regions
	# Output: these parameters in a hashmap format

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	# params
	config = {}
	for row in reader:
		if len(row)==0: continue
		if 'folder' in row[0]: break
		if '#' in row[0] or row[0]=='': continue
		config[row[0]] = row[1]

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

	return config, folders, regions

def processFolders(parentDir, folders, regions, numChr):
	folderToGraph = {}
	for folder in folders:

		# setting up folders
		os.chdir(parentDir)
		if not glob.glob(folder + '/'): os.mkdir(folder)
		os.chdir(folder)
		
		if not glob.glob('bins/'): os.system("mkdir bins")
		os.chdir('bins')
		#os.system("rmdir *") # removing empty directories
		regionFolders = ' '.join(regions.keys())
		try: subprocess.check_output("mkdir " + regionFolders, stderr = a, shell=True)
		except: pass

		# splitting up bedgraph if not done already
		os.chdir("..")
		if not glob.glob("bedGraphByChr/"): os.system("mkdir bedGraphByChr")
		
		os.chdir("bedGraphByChr")
		if len(glob.glob("*.bedGraph")) != numChr + 2: # 23 --> 25 for human (chr1-22, x, y, m), 20 --> 22 for mouse (chr1-19, x, y, m)
			os.system("rm -f *.bedGraph")
			print "Splitting up bedgraph for " + folder
			#SNF mod: awk -> gawk 
			cmd = "gawk '{print >> $1\".bedGraph\"}' " + folders[folder][0]
			print cmd
			os.system(cmd)

		# making folder to bedgraph relationship
		binFolder = parentDir + '/' + folder + '/bins/'
		graphFolder = parentDir + '/' + folder + '/bedGraphByChr/'
		folderToGraph[folder] = [binFolder, graphFolder, folders[folder][1]] #bin folder --> [graph folder, strand]

	return folderToGraph

# checks whether file is a bed file
def isBed(row):
	if 'chr' not in row[0]: return False
	try:
		a=int(row[1])
		b=int(row[2])
	except: return False
	return True
	
def getChrToRegion(fn, header):
	with open(fn, 'r') as ifile:	
		regions = collections.defaultdict(lambda: []) # by chromosome
		if header: ifile.next()
		
		counter = 0
		for line in ifile:
			counter += 1
			row = line.split() # added by SNF 
			if not isBed(row):
				print "Line %d of file %s is not in BED format" %(counter, fn)
				exit()
			regions[row[0]].append(row)
	return regions

def processRegions(regions):
	regionToChrMap = {}
	for region in regions:
		info = regions[region]
		loc = info[0]
		isHeader = True if info[1] == 'y' else False
		regionToChrMap[region] = getChrToRegion(loc, isHeader)

	return regionToChrMap

def main():
	# check
	if len(sys.argv) < 2: 
		print "Need configuration file."
		exit()

	# log file
	# logfile = open('logs/' + str(random.randrange(1,1000)) + '.log', 'w')
	# logwriter = csv.writer(logfile, 'textdialect')
	
	# reading config file
	config, folders, regions = readConfigFile(sys.argv[1])
	print "\nRead configuration file"
	for c in config: print c + ':', config[c]

	# chromosome configuration
	organism = config['organism(mm9 or hg19)']
	numChr = 23 if organism == 'hg19' else 20
	allChroms = ['chr' + str(x) for x in range(1,numChr)]
	allChroms.extend(['chrX', 'chrY', 'chrM'])
	threads = int(config['threads'])
	numProcs = threads
	
	# processing folders and bedgraphs
	parentDir = config["parentDir"]
	folderToGraph = processFolders(parentDir, folders, regions, numChr)
	print "\nProcessed folders: ", ', '.join(folderToGraph.keys())
	for f in folderToGraph:
		print f + ':', folderToGraph[f][0]
	
	# processing regions, checking that they are valid bed files
	regionToChrMap = processRegions(regions)
	print "\nProcessed regions: ", ', '.join(regionToChrMap.keys())
	for r in regionToChrMap:
		print r, regions[r][0]

	# making bins
	print "\nReading in bedgraphs and making profiles for each region"
	binLength = int(sys.argv[2]) # how long is the bin where we put bedgraph regions?
	#binLength = 2000000
	xstart = datetime.now()
	for folder in folderToGraph:
		# if my bedgraph is stranded and my regions are stranded, only 
		# use the regions that correspond to the bedgraph strand
		[binFolder, graphFolder, folderStrand] = folderToGraph[folder]

		# process all regions for each sub-bedgraph
		for i in range(len(allChroms)):
			chroms = allChroms[(numProcs*i):(numProcs*(i+1))]			
			reads = readBedGraph(graphFolder, chroms, binLength)
			if reads == {}: continue
			
			tstart = datetime.now()
			for region in regions:
				info = regions[region]
				stranded = True if info[2]=='y' else False
				limitSize = True if info[3]=='y' else False
				extendRegion = True if info[5]=='y' else False
				numBins = int(info[4])
				
				regionProcess(binFolder, region, regionToChrMap[region], chroms, stranded, folderStrand, limitSize, numBins, extendRegion, reads, binLength)
			tend = datetime.now()
			delta = tend - tstart
			print delta.total_seconds()
			
	xend = datetime.now()
	delta = xend - xstart
	print delta.total_seconds()

	# merging bins for each chromosome, then make metagene
	folders = folderToGraph.keys() 
	numPerProc = len(folders)/numProcs + 1 if len(folders) > numProcs else 1
	numJobs = numProcs if (numProcs < len(folders)) else len(folders)
	procs = []

	for i in range(numJobs): 
		p = multiprocessing.Process(target=folderWorker, args=(i * numPerProc, (i + 1) * numPerProc, folders, folderToGraph, regions))
		procs.append(p)
		p.start()
	for p in procs: p.join()
	print "Made metagenes"

	# merging all files, and writing average files
	name = config["name"]
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	os.chdir(parentDir)
	os.system("mkdir averages")
	for region in regions:
		for folder in folderToGraph:
			binFolder = folderToGraph[folder][0]
			os.chdir(binFolder + '/' + region + '/')
			fn = "avgraw_" + folder + "_" + region 
			regionToFolderAvgs[region][folder] = processFile(fn)
		writeFile(name + '_' + region, regionToFolderAvgs[region], parentDir + '/averages/')

if __name__ == '__main__':
	main()

