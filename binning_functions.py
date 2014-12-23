# binning_functions.py
# 8/29/14
# helper functions to make bins for each region

import os, glob, csv, re, collections, math, multiprocessing, sys, time, logging
from datetime import datetime

logger = logging.getLogger('')

# gets the average score within this bin
# end of bin is NOT inclusive, because end of read in bedgraph is also NOT inclusive
def getAverageScore(readsForChrom, currStart, currEnd, binNum, rn):
	totalLength = currEnd - currStart
	totalScore = 0
	readNumber = rn 			# which read to start at

	while True:		
		try: 
			readList = readsForChrom[binNum]
			read = readList[readNumber]
		except: #reached end of bin
			binNum += 1
			readNumber = 0
			if binNum not in readsForChrom: break  #reached end of chromosome
			readList = readsForChrom[binNum]
			read = readList[readNumber]
			
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
			totalScore += score * (currEnd - currStart)
			#print 'totalScore', totalScore, currEnd - currStart + 1
			break				# don't advance read number because the read extends past end of gene
		else: # end of gene is past the read
			totalScore += score * (readEnd - currStart)
			#print 'totalScore', totalScore, readEnd - currStart + 1
			readNumber += 1 	# advance read number because this read won't be needed anymore

		currStart = readEnd + 1

	return float(totalScore)/totalLength, readNumber, binNum

def getInitialReadNumber(readsForChrom, binNum, currStart):
	# better way to find the initial readnumber; saves a few seconds
	readNumber = 0 #index of read, updated each time
	readList = readsForChrom[binNum]
	left = 0
	right = len(readList) - 1	
	curr = 0
	while True:
		curr = (left + right)/2
		#print curr, currStart, readList[curr]
		#time.sleep(0.3)
		if curr == left or curr == right: break
		if currStart < readList[curr][0]: right = curr
		else: left = curr	

	return curr


# gets the array of bins for each gene
# end of gene is NOT inclusive
def getBins(start, end, numBins, readsForChrom, binLength):

	# walking through each bin for the region
	scores = []
	readNumber = getInitialReadNumber(readsForChrom, binNum, currStart) - 1
	currStart = start
	binNum = start/binLength    	
	spacingPerBin = int(math.ceil((end - start)/float(numBins)))

	while currStart < end:
		currEnd = currStart + spacingPerBin # end of my window
		if currEnd > end: currEnd = end # last bin can't go past TES
		
		score, readNumber, binNum = getAverageScore(readsForChrom, currStart, currEnd, binNum, readNumber) #updates read number also
		scores.append(score)
		currStart = currEnd # new start of my window

	return scores

# for each chromosome, get bins corresponding to each region in the chromosome
def regionWorker(binFolder, regionType, chrom, chrToIndivRegions, stranded, folderStrand, limitSize, numBins, extendRegion, readsForChrom, binLength):
	ofile = open(binFolder + "/" + regionType + "/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for region in chrToIndivRegions[chrom]:
		start = int(region[1])
		end = int(region[2])
		
		#print "SNF: regionworker start %d end %d limitSize %r region %s" % (start, end, limitSize, region)

		# binning doesn't really make sense here so just ignore
		if limitSize:
			if end - start < 200 or end - start > 200000: continue 

		# extending to 1x upstream and downstream
		if extendRegion:
			length = end - start
			start = start - length
			end = end + length
			region[1] = start
			region[2] = end
			
		# only taking the regions that match the strand of bedgraph, 
		# if bedgraph and region file are both stranded
		if stranded: strand = region[5]
		#if folderStrand != '0' and stranded and strand != folderStrand: continue 

		# getting bins and reading from end to start if region is antisense
		regionBins = getBins(start, end, numBins, readsForChrom, binLength)
		if stranded and strand == '-': regionBins = regionBins[::-1] 

		outputRow = region
		outputRow.append(sum(regionBins)) 
		outputRow.extend(regionBins)
		writer.writerow(outputRow)

	ofile.close()
	

# gets the array of bins for each gene
# start and end of each block are NOT inclusive
def blockGetBins(blocks, numBins, readsForChrom, binLength):
		
	# variables that will not change
	scores = []
	totalLength = sum([bk[1]-bk[0] for bk in blocks])
	end = blocks[-1][1]
	spacingPerBin = int(math.ceil(totalLength/float(numBins)))

	# blocks is [[block1s, block1e], [block2s, block2e], ...]
	# variables that will change as the loop iterates
	binNum = currStart/binLength    				
	currStart = blocks[0][0]
	readNumber = getInitialReadNumber(readsForChrom, binNum, currStart) - 1
		
	blockNum = 0
	contBlock = False
	oldPartialLength = 0
	oldPartialScore = 0

	while True:
		blockStart, blockEnd = blocks[blockNum][0], blocks[blockNum][1] #info about current block
		currEnd = currStart + spacingLeft
			
		if currEnd >= end: currEnd = end # last bin can't go past TES
		if currEnd <= blockEnd:
			print "case 1: ", currStart, currEnd
			score,readNumber,binNum = getAverageScore(readsForChrom, currStart, currEnd, binNum, readNumber)
			newPartialLength = currEnd - currStart
			score = (score * newPartialLength + partialScore * oldPartialLength)/float(spacingPerBin) #weighted average
			scores.append(score)
			
			# reset things
			contBlock = False
			spacingLeft = spacingPerBin
			oldPartialLength = 0
			oldPartialScore = 0
			
			# define start for next iteration of while loop
			if currEnd < blockEnd: currStart = currEnd # new start of my window
			else: 
				if currEnd == end: break
				blockNum += 1
				currStart = blocks[blockNum][0]
		else:
			print "case 2: ", currStart, blockEnd, spacingLeft\
			currEnd = blockEnd
			score,readNumber,binNum = getAverageScore(readsForChrom, currStart, currEnd, binNum, readNumber)
			newPartialLength = currEnd - currStart
			score = (score * newPartialLength + partialScore * oldPartialLength)/float(newPartialLength + oldPartialLength) #weighted average
			
			# carry over to new block
			contBlock = True
			oldPartialScore = score
			oldPartialLength = oldPartialLength + newPartialLength		# remaining bin length is reduced
			blockNum += 1
			
			# define start for next iteration of while loop
			currStart = blocks[blockNum][0] # start at the next bin
		
	return scores

# for each chromosome, get bins corresponding to each region in the chromosome
def blockRegionWorker(binFolder, regionType, chrom, chrToIndivRegions, stranded, folderStrand, limitSize, numBins, extendRegion, readsForChrom, binLength):
	ofile = open(binFolder + "/" + regionType + "/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for region in chrToIndivRegions[chrom]:
	
		# get info for each region
		start = int(region[1])
		end = int(region[2])
		blockLengthStr = region[10] #fix
		blockStartStr = region[11] #fix
		blockLengths = [int(x) for x in blockLengthStr.split(',')]
		blockStarts = [int(x) for x in blockStartStr.split(',')]
		
		# make blocks
		blocks = []
		for i in range(len(blockStarts)):
			blockStart = start + blockStarts[i]
			blockEnd = blockStart + blockLengths[i]
			blocks.append([blockStart, blockEnd])
		
		# limits
		totalLength = sum(blockLengths)
		if limitSize:
			if totalLength < 200 or totalLength > 200000: continue 
		
		# getting bins and reading from end to start if region is antisense
		# only taking the regions that match the strand of bedgraph, 
		# if bedgraph and region file are both stranded
		regionBins = blockGetBins(blocks, numBins, readsForChrom, binLength)
		if stranded: strand = region[5]
		if stranded and strand == '-': regionBins = regionBins[::-1] 

		# write to file
		outputRow = region
		outputRow.append(sum(regionBins)) 
		outputRow.extend(regionBins)
		writer.writerow(outputRow)

	ofile.close()
	
#Loading reads for each bedgraph, for each chromosome
def getReads(chrom, graph, binLength):
	ifile = open(graph, 'r')
	reader = csv.reader(ifile, 'textdialect')
	readsForChrom = {}

	for row in reader:
		start = int(row[1])
		end = int(row[2])
		score = float(row[3])

		# which bins does each interval go into?
		binNum = start/binLength
		bin1 = binNum
		bin2 = bin1 - 1
		
		# result: bins are overlapping
		if bin1 not in readsForChrom.keys(): readsForChrom[bin1] = []
		readsForChrom[bin1].append([start, end, score])
		if bin2 not in readsForChrom.keys(): readsForChrom[bin2] = []
		readsForChrom[bin2].append([start, end, score])

	ifile.close()
	return readsForChrom

# reads bedgraph and does region processing for each chromosome
def processEachChrom(chrom, binFolder, graphFolder, folderStrand, binLength, regions, regionToChrMap):
	# reads in bedgraph
	if not glob.glob(graphFolder + chrom + '.bedGraph'): return
	logger.info('%s %s', chrom, graphFolder + chrom + '.bedGraph')
	readsForChrom = getReads(chrom, graphFolder + chrom + '.bedGraph', binLength)
	
	# processes regions
	# tstart = datetime.now()
	for region in regions:
		info = regions[region]
		stranded = True if info[2]=='y' else False
		limitSize = True if info[3]=='y' else False
		extendRegion = True if info[5]=='y' else False
		numBins = int(info[4])
		
		chrToIndivRegions = regionToChrMap[region]
		regionWorker(binFolder, region, chrom, chrToIndivRegions, stranded, folderStrand, limitSize, numBins, extendRegion, readsForChrom, binLength)
	
	logger.info('%s done: %s', chrom, ', '.join(regions.keys()))
	# tend = datetime.now()
	# delta = tend - tstart
	# print delta.total_seconds()