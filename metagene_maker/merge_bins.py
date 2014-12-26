# merge_bins.py
# 8/29/14
# merges bins into metagene for each region type for each bedgraph

import os, glob, csv, re, multiprocessing, logging

logger = logging.getLogger('')

# concats chr1, chr2, etc --> allchr_sorted.txt
def concatChrs(start, end, folders, folderToGraph, regions):
	for i in range(start, end):
		try: folder = folders[i]
		except: continue
		binFolder = folderToGraph[folder][0]
		for region in regions: 
			folder = binFolder + '/' + region + '/'
			os.chdir(folder)
			if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")
			os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")
			os.system("rm -f allchr.txt")

# runs RScript on allchr_sorted.txt
def runRScript(start, end, folders, folderToGraph, regions):
	script = os.path.dirname(__file__) + "/rscript/makeMetagenePlot.r"
	
	for i in range(start, end):
		try: folder = folders[i]
		except: continue
		binFolder = folderToGraph[folder][0]
		for region in regions: 
			os.chdir(binFolder + '/' + region + '/')

			info = regions[region]
			numBins = info[2]
			rcmd = ' '.join(["Rscript", script, folder, region, "6", "3", numBins, "allchr_sorted.txt"])
			logger.info('%s', rcmd)
			os.system(rcmd)

# distill (+) and (-) into sense and antisense
def processPaired(pair, folderPairs, regions, folderToGraph, parentDir):
	plusSample = folderPairs[pair][0]
	plusParentDir = parentDir + '/' + plusSample + '/bins/'
	minusSample = folderPairs[pair][1]
	minusParentDir = parentDir + '/' + minusSample + '/bins/'
	
	for region in regions: #we are guaranteed that the regions must be stranded, so they are either + or -
		logger.info("\nSplitting %s %s", pair, region)
		# split allchr into plus and minus
		plusdir = plusParentDir + region + '/'
		os.chdir(plusdir)
		os.system("gawk -F '\t' '{print >> \"allchr_\" $6 \".txt\"}' allchr_sorted.txt")
		minusdir = minusParentDir + region + '/'
		os.chdir(minusdir)
		os.system("gawk -F '\t' '{print >> \"allchr_\" $6 \".txt\"}' allchr_sorted.txt")
			
		# combine plus plus, minus minus for sense
		logger.info("Making sense for %s %s", pair, region) 
		opath = parentDir + "/" + pair + "_sense/bins/"
		os.chdir(opath)
		if not glob.glob(region): os.system("mkdir " + region)
		os.chdir(region)
		file1 = plusdir + "/allchr_+.txt"
		file2 = minusdir + "/allchr_-.txt"
		os.system("cat " + file1 + " " + file2 + " > " + "allchr_sorted.txt")
		
		# combine plus minus, minus plus for antisense
		logger.info("Making antisense for %s %s", pair, region)
		opath = parentDir + "/" + pair + "_antisense/bins/"
		os.chdir(opath)
		if not glob.glob(region): os.system("mkdir " + region)
		os.chdir(region)
		file1 = minusdir + "/allchr_+.txt"
		file2 = plusdir + "/allchr_-.txt"
		os.system("cat " + file1 + " " + file2 + " > " + "allchr_sorted.txt")
		
		
		
		
		
			
# take in an avg_*_* file and extract the bin values
def processFile(fileName, isMinus):
	avgFile = glob.glob(fileName + "*.txt")[0]
	logger.info(avgFile)
	ifile = open(avgFile, 'r')
	reader = csv.reader(ifile, 'textdialect')

	values = []
	reader.next()
	if isMinus:
		for row in reader:
			values.append(-float(row[1]))
	else:
		for row in reader:
			values.append(row[1])

	ifile.close()
	return values

# write the mapping to a file
def writeFile(name, mapping, direc):
	ofile = open(direc + name + '.txt', 'w')
	writer = csv.writer(ofile, 'textdialect')
	for treatment in mapping:
		outputRow = [treatment]
		outputRow.extend(mapping[treatment])
		writer.writerow(outputRow)

	ofile.close()
