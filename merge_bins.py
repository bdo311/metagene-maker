# merge_bins.py
# 8/29/14
# merges bins into metagene for each region type for each bedgraph

import os, glob, csv, re, multiprocessing
csv.register_dialect("textdialect", delimiter='\t')

def mergeChrAndSort(folder):
	os.chdir(folder)
	if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")
	os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")

def mergeChr(folder):
	os.chdir(folder)
	if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")

script = "makeMetagenePlot.r"
def folderWorker(start, end, folders, folderToGraph, regions):
	folder = folders[i]
	binFolder = folderToGraph[folder][0]
	for region in regions: 
		os.chdir(binFolder + '/' + region + '/')
		mergeChr(binFolder + '/' + region + '/')

		info = regions[region]
		numCols, nameCol, numBins = int(row[8]), int(row[3]), int(info[10])
		os.system(' '.join(["Rscript", script, folder, region, numCols, nameCol, numBins]))

# take in an avg_*_* file and extract the bin values
def processFile(fileName):
	avgFile = glob.glob(fileName + "*.txt")[0]
	print avgFile
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
