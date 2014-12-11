# binning_functions.py
# 8/29/14
# helper functions to make bins for each region

import os, glob, csv, re, collections, math, multiprocessing, sys

# gets the average score within this bin
def getAverageScore(chrom, currStart, currEnd, binNum, rn, reads):
	totalLength = currEnd - currStart + 1
	totalScore = 0
	readNumber = rn 			# which read to start at

	if binNum not in reads[chrom]: return 0, readNumber
	readList = reads[chrom][binNum][rn:]
	for read in readList:
		readStart, readEnd, score = read[0], read[1], read[2]
		#print 'read', readStart, readEnd, score

		if readEnd < currStart:
			readNumber += 1 	# never read that read again
			continue 			# moves on to next read
		if currEnd < readStart: 
			#print 'totalScore', totalScore
			break 				# end of gene < start of read. No coverage, so score is 0

		if currStart < readStart: currStart = readStart # ignore zeros that are in the read region
		if readEnd >= currEnd: 	# end of gene is within the read
			totalScore += score * (currEnd - currStart + 1)
			#print 'totalScore', totalScore, currEnd - currStart + 1
			break				# don't advance read number because the read extends past end of gene
		else: # end of gene is past the read
			totalScore += score * (readEnd - currStart + 1)
			#print 'totalScore', totalScore, readEnd - currStart + 1
			readNumber += 1 	# advance read number because this read won't be needed anymore

		currStart = readEnd + 1

	return float(totalScore)/totalLength, readNumber

# gets the array of bins for each gene
def getBins(chrom, start, end, numBins, rn, reads):
	spacingPerBin = int(math.ceil((end - start)/float(numBins)))
	binNum = start/500000    

	scores = []
	currStart = start
	readNumber = 0 #index of read, updated each time
	while currStart < end:
		currEnd = currStart + spacingPerBin - 1 # end of my window
		if currEnd > end: currEnd = end # last bin can't go past TES
		
		score, readNumber = getAverageScore(chrom, currStart, currEnd, binNum, readNumber, reads) #updates read number also
		scores.append(score)
		currStart = currEnd + 1 # new start of my window

	return scores

# for each chromosome, get bins corresponding to each region in the chromosome
def regionWorker(binFolder, regionType, chrom, chrToRegion, startCol, endCol, stranded, folderStrand, strandCol, limitSize, numBins, extendRegion, reads):
	ofile = open(binFolder + "/" + regionType + "/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for region in chrToRegion[chrom]:
		start = int(region[startCol])
		end = int(region[endCol])
		
		#print "SNF: regionworker start %d end %d limitSize %r region %s" % (start, end, limitSize, region)

		# binning doesn't really make sense here so just ignore
		if limitSize:
			if end - start < 200 or end - start > 200000: continue 

		# extending to 1x upstream and downstream
		if extendRegion:
			length = end - start
			start = start - length
			end = end + length
			region[startCol] = start
			region[endCol] = end
		# only taking the regions that match the strand of bedgraph, 
		# if bedgraph and region file are both stranded
		strand = region[strandCol]
		#if folderStrand != '0' and stranded and strand != folderStrand: continue 

		# getting bins and reading from end to start if region is antisense
		regionBins = getBins(chrom, start, end, numBins, 0, reads)
		if stranded and strand == '-': regionBins = regionBins[::-1] 

		outputRow = region
		outputRow.append(sum(regionBins)) 
		outputRow.extend(regionBins)
		writer.writerow(outputRow)

	ofile.close()
	
def regionProcess(binFolder, regionType, chrToRegion, chroms, startCol, endCol, stranded, folderStrand, strandCol, limitSize, numBins, extendRegion, reads):
	print "Working on " + regionType
	procs = []
	for chrom in chroms:
		p = multiprocessing.Process(target=regionWorker, args=(binFolder, regionType, chrom, chrToRegion, startCol, endCol, stranded, folderStrand, strandCol, limitSize, numBins, extendRegion, reads))
		procs.append(p)
		p.start()
	for p in procs:
		p.join()	# wait till all are done before moving on
		











		
# Loading reads for each bedgraph, for each chromosome
def getReads(readQueue, chrom, graph):
	ifile = open(graph, 'r')
	reader = csv.reader(ifile, 'textdialect')
	readsForChrom = {}

	for row in reader:
		start = int(row[1])
		end = int(row[2])
		score = float(row[3])

		# which bins does each interval go into?
		multipleOf500k = start/500000
		bin1 = multipleOf500k
		bin2 = bin1 - 1
		
		# result: bins are overlapping
		if bin1 not in readsForChrom.keys(): readsForChrom[bin1] = []
		readsForChrom[bin1].append([start, end, score])
		if bin2 not in readsForChrom.keys(): readsForChrom[bin2] = []
		readsForChrom[bin2].append([start, end, score])

	ifile.close()
	readQueue.put((chrom, readsForChrom))
	print chrom, 'done'

def readBedGraph(ifolder, chroms):
	manager = multiprocessing.Manager()
	readQueue = manager.Queue()

	# get all chromosomes in a queue
	procs = []
	for chrom in chroms:
		print chrom, ifolder + chrom + '.bedGraph'
		p = multiprocessing.Process(target=getReads, args=(readQueue, chrom, ifolder + chrom + '.bedGraph'))
		p.start()
		procs.append(p)
	for proc in procs: proc.join()

	# make a large dictionary
	print "Making dictionary"
	reads = {}
	while not readQueue.empty():
		readTuple = readQueue.get()
		reads[readTuple[0]] = readTuple[1]
	return reads
