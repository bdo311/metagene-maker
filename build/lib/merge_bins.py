# merge_bins.py
# 8/29/14
# merges bins into metagene for each region type for each bedgraph

import os, glob, csv, re, multiprocessing, logging

logger = logging.getLogger('')

def mergeChrAndSort(folder):
	os.chdir(folder)
	if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")
	os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")

def mergeChr(folder):
	os.chdir(folder)
	#logger.info(folder)	
	if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")

def folderWorker(start, end, folders, folderToGraph, regions):
	script = os.path.join(os.path.dirname(__file__), "/makeMetagenePlot.r")
	
	for i in range(start, end):
		try: folder = folders[i]
		except: continue
		binFolder = folderToGraph[folder][0]
		for region in regions: 
			os.chdir(binFolder + '/' + region + '/')
			mergeChr(binFolder + '/' + region + '/')

			info = regions[region]
			numBins = info[4]
			rcmd = ' '.join(["Rscript", script, folder, region, "6", "3", numBins, "allchr.txt"])
			logger.info('%s', rcmd)
			os.system(rcmd)

# take in an avg_*_* file and extract the bin values
def processFile(fileName):
	avgFile = glob.glob(fileName + "*.txt")[0]
	logger.info(avgFile)
	ifile = open(avgFile, 'r')
	reader = csv.reader(ifile, 'textdialect')

	values = []
	reader.next()
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
