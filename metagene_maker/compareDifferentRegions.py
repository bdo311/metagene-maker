# compareDifferentRegions.py
# get lists of average read densities for different regions across different samples based on a config file

import csv, sys, os, glob, multiprocessing
import numpy as np
csv.register_dialect("textdialect",delimiter='\t')


def readConfig(fn):
	samples = []
	regions = {} # regions --> [start, end] inclusive, 1-based

	with open(fn,'r') as ifile:
		reader = csv.reader(ifile, 'textdialect')
		for row in reader:
			if 'region' in row[0]: break
			if len(row)==0 or '#' in row[0]: continue
			name = row[0]
			if name.rstrip() == '': continue
			if not glob.glob(name + '/'):
				print "Sample " + name + " not found in folder. Exiting."
				exit()
			samples.append(name)
			
		for row in reader:
			if len(row)==0 or '#' in row[0]: continue
			if row[0]=='': continue
			try:
				print row
				a,b=int(row[1]),int(row[2])
				regions[row[0]] = [a,b]
			except:
				print row[0] + " has invalid configuration. Exiting."
	
	return samples, regions		
		
def getSortedValues(sample, region, config):
	[start, end] = config
	fn = sample + "/bins/" + region + '/allchr_sorted.txt'

	if not glob.glob(fn):
		print fn + " not found. Exiting."
		exit()
		
	nameToValues = {}
	with open(fn, 'r') as ifile:
		reader = csv.reader(ifile, 'textdialect')
		counter = 0
		for row in reader:
			counter += 1
			#if counter == 10: break
			name = '__'.join(row[:6])
			values = [float(x) for x in row[(7+start-1):(7+end)]]
			avg = np.mean(values)
			nameToValues[name] = avg
		
	keys = nameToValues.keys()
	keys.sort()
	sortedValues = []
	for key in keys:
		sortedValues.append(nameToValues[key])
	
	return keys, sortedValues	
	
	
def regionWorker(start, end, regionNames, regions, samples, suffix):
	for region in regionNames[start:end]:
		print "Processing " + region
		regionToSamples = {}
		keys = []
		for sample in samples:
			keys, regionToSamples[sample] = getSortedValues(sample, region, regions[region])
			
		
		ofile = open("regionLists/" + region + "_list_" + suffix + ".txt",'w')
		writer = csv.writer(ofile, 'textdialect')
		header = ['regionName']
		header.extend(sorted(regionToSamples.keys()))
		writer.writerow(header)
		numRegions = len(regionToSamples[header[1]])
		for i in range(numRegions):
			outputRow = [keys[i]]
			outputRow.extend([regionToSamples[x][i] for x in header[1:]])
			writer.writerow(outputRow)
			
		ofile.close()
		
def main():
	if len(sys.argv) != 3:
		print "Usage: compareDifferentRegions <config file> <suffix>"
		exit()
		
	configFile = sys.argv[1]
	suffix = sys.argv[2]
	samples, regions = readConfig(configFile)
	if not glob.glob("regionLists/"): os.system("mkdir regionLists")
	
	numProcs = 4
	numPerProc = len(regions)/numProcs + 1 if len(regions) > numProcs else 1
	numJobs = numProcs if (numProcs < len(regions)) else len(regions)
	procs = []
	
	regionNames = regions.keys()
	for i in range(numJobs): 
		p = multiprocessing.Process(target=regionWorker, args=(i * numPerProc, (i + 1) * numPerProc, regionNames, regions, samples, suffix))
		procs.append(p)
		p.start()
	for p in procs: p.join()

if __name__=='__main__':
	main()
