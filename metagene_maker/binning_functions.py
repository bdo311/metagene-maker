# binning_functions.py
# 8/29/14, last updated 12/28/14
# helper functions to make bins for each region

import os, glob, csv, re, collections, math, multiprocessing, sys, time, logging
import numpy as np
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
			if readsForChrom[binNum] == []: break  #reached end of chromosome
			readList = readsForChrom[binNum]
			read = readList[readNumber]
			
		#print read
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

	#print totalScore
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
		if curr == left or curr == right: break
		if currStart < readList[curr][0]: right = curr
		else: left = curr	

	return curr


# gets the array of bins for each gene
# end of gene is NOT inclusive
def getBins(start, end, numBins, readsForChrom, binLength):
	# walking through each bin for the region
	scores = []
	binNum = start/binLength   
	currStart = start
	readNumber = getInitialReadNumber(readsForChrom, binNum, currStart)
	if readNumber < 0: readNumber = 0
	spacingPerBin = int(math.ceil((end - start)/float(numBins)))

	while currStart < end:
		currEnd = currStart + spacingPerBin # end of my window
		if currEnd > end: currEnd = end # last bin can't go past TES
		
		#print currStart, currEnd, binNum, readNumber
		score, readNumber, binNum = getAverageScore(readsForChrom, currStart, currEnd, binNum, readNumber) #updates read number also
		scores.append(score)
		currStart = currEnd # new start of my window
	
	# make sure the length of scores is right
	if len(scores) == 0: scores = [0] * numBins #if start > end (when calculating TR for genes < 300bp, scores will be empty, and we need to fix that or else interp won't work
	if len(scores) != numBins:
		a=map(lambda x: float(x)*len(scores)/numBins, range(numBins)) # convert my desired scale to the current scale
		
		scores = np.interp(a, range(len(scores)), scores)
		
	assert len(scores)==numBins
	return list(scores)

# for each chromosome, get bins corresponding to each region in the chromosome
def regionWorker(binFolder, regionType, chrom, chrToIndivRegions, limitSize, numBins, extension, sideNumBins, readsForChrom, binLength):
	ofile = open(binFolder + "/" + regionType + "/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for region in chrToIndivRegions[chrom]:
		start = int(region[1])
		end = int(region[2])
		
		# binning doesn't really make sense here so just ignore
		if limitSize:
			if end - start < 200 or end - start > 200000: continue 

		# getting bins and reading from end to start if region is antisense
		if extension > 0:
			leftSideBins = getBins(start - extension, start, sideNumBins, readsForChrom, binLength)
			regionBins = getBins(start, end, numBins, readsForChrom, binLength)
			rightSideBins = getBins(end, end + extension, sideNumBins, readsForChrom, binLength)
			
			leftSideBins.extend(regionBins)
			leftSideBins.extend(rightSideBins)
			regionBins = leftSideBins
		else: 
			regionBins = getBins(start, end, numBins, readsForChrom, binLength)
			
	
		strand = region[5]
		if strand == '-': regionBins = regionBins[::-1] 

		outputRow = region[:6]
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
	currStart = blocks[0][0]
	binNum = currStart/binLength    				
	readNumber = getInitialReadNumber(readsForChrom, binNum, currStart) - 1
	if readNumber < 0: readNumber = 0
		
	blockNum = 0
	contBlock = False
	oldPartialLength = 0
	oldPartialScore = 0

	while True:
		blockStart, blockEnd = blocks[blockNum][0], blocks[blockNum][1] #info about current block
		currEnd = currStart + spacingPerBin - oldPartialLength
			
		if currEnd >= end: currEnd = end # last bin can't go past TES
		if currEnd <= blockEnd:
			# print "case 1: ", currStart, currEnd, oldPartialLength, contBlock
			score,readNumber,binNum = getAverageScore(readsForChrom, currStart, currEnd, binNum, readNumber)
			newPartialLength = currEnd - currStart
			# print 'individual score: ', score
			# print 'newPartialLength: ', newPartialLength
			# print 'oldPartialScore: ', oldPartialScore
			# print 'oldPartialLength: ', oldPartialLength
			score = (score * newPartialLength + oldPartialScore * oldPartialLength)/float(newPartialLength + oldPartialLength) #weighted average
			# print 'overall score: ', score
			scores.append(score)
			
			# reset things
			contBlock = False
			oldPartialLength = 0
			oldPartialScore = 0
			
			# define start for next iteration of while loop
			if currEnd < blockEnd: currStart = currEnd # new start of my window
			else: 
				if currEnd == end: break
				blockNum += 1
				currStart = blocks[blockNum][0]
		else:
			# print "case 2: ", currStart, blockEnd, oldPartialLength, contBlock
			currEnd = blockEnd
			score,readNumber,binNum = getAverageScore(readsForChrom, currStart, currEnd, binNum, readNumber)
			newPartialLength = currEnd - currStart
			# print 'individual score: ', score
			# print 'newPartialLength: ', newPartialLength
			# print 'oldPartialScore: ', oldPartialScore
			# print 'oldPartialLength: ', oldPartialLength
			
			score = (score * newPartialLength + oldPartialScore * oldPartialLength)/float(newPartialLength + oldPartialLength) #weighted average
			# print 'overall score: ', score

			# carry over to new block
			contBlock = True
			oldPartialScore = score
			oldPartialLength = oldPartialLength + newPartialLength		# remaining bin length is reduced
			blockNum += 1
			
			# define start for next iteration of while loop
			currStart = blocks[blockNum][0] # start at the next bin
		
	# make sure the length of scores is right
	if len(scores) == 0: scores = [0] * numBins #if start > end (when calculating TR for genes < 300bp, scores will be empty, and we need to fix that or else interp won't work
	if len(scores) != numBins:
		a=map(lambda x: float(x)*len(scores)/numBins, range(numBins)) # convert my desired scale to the current scale
		scores = np.interp(a, range(len(scores)), scores)
	assert len(scores)==numBins
	return scores

# for each chromosome, get bins corresponding to each region in the chromosome
def blockRegionWorker(binFolder, regionType, chrom, chrToIndivRegions, limitSize, numBins, extension, sideNumBins, readsForChrom, binLength):
	ofile = open(binFolder + "/" + regionType + "/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for region in chrToIndivRegions[chrom]:
	
		# get info for each region
		start = int(region[1])
		end = int(region[2])
		blockLengthStr = region[10] #fix
		blockStartStr = region[11] #fix
		
		# split block info
		blockLengths = []
		for x in blockLengthStr.split(','):
			try: blockLengths.append(int(x))
			except: pass
		blockStarts = []
		for x in blockStartStr.split(','):
			try: blockStarts.append(int(x))
			except: pass
				
		# make blocks
		blocks = []
		for i in range(len(blockStarts)):
			blockStart = start + blockStarts[i]
			blockEnd = blockStart + blockLengths[i]
			blocks.append([blockStart, blockEnd])
		# print blocks
		
		# limits
		totalLength = sum(blockLengths)
		if limitSize:
			if totalLength < 200 or totalLength > 200000: continue 
			
		# getting bins and reading from end to start if region is antisense
		# only taking the regions that match the strand of bedgraph, 
		# if bedgraph and region file are both stranded
		if extension > 0:
			leftSideBins = getBins(start - extension, start, sideNumBins, readsForChrom, binLength)
			regionBins = blockGetBins(blocks, numBins, readsForChrom, binLength)
			rightSideBins = getBins(end, end + extension, sideNumBins, readsForChrom, binLength)
			
			leftSideBins.extend(regionBins)
			leftSideBins.extend(rightSideBins)
			regionBins = leftSideBins
		else: 
			regionBins = blockGetBins(blocks, numBins, readsForChrom, binLength)
			
		strand = region[5]
		if strand == '-': regionBins = regionBins[::-1] 

		# write to file
		outputRow = region[:6]
		outputRow.append(sum(regionBins)) 
		outputRow.extend(regionBins)
		writer.writerow(outputRow)

	ofile.close()
	
#Loading reads for each bedgraph, for each chromosome
def getReads(chrom, graph, binLength):
	ifile = open(graph, 'r')
	reader = csv.reader(ifile, 'textdialect')
	readsForChrom = collections.defaultdict(lambda: [])

	for row in reader:		
		start, end, score = int(row[1]), int(row[2]), float(row[3])

		# which bins does each interval go into?
		bin1 = start/binLength
		readsForChrom[bin1].append([start, end, score])

		bin2 = end/binLength
		while bin1 < bin2:
			bin1 = bin1 + 1
			readsForChrom[bin2].append([start, end, score])
			
	ifile.close()
	return readsForChrom

# reads bedgraph and does region processing for each chromosome
def processEachChrom(chrom, binFolder, graphFolder, binLength, regions, regionToChrMap, regionToBedType):
	# reads in bedgraph
	if not glob.glob(graphFolder + chrom + '.bedGraph'): return

	logger.info('%s %s', chrom, graphFolder + chrom + '.bedGraph')
	readsForChrom = getReads(chrom, graphFolder + chrom + '.bedGraph', binLength)
	#print readsForChrom
	#logger.info('Read %s', chrom)
	# processes regions
	for region in regions:
		info = regions[region]
		limitSize = True if info[1]=='y' else False
		sideExtension = int(info[3])
		sideNumBins = int(info[4])
		numBins = int(info[2])
		
		chrToIndivRegions = regionToChrMap[region]
		if regionToBedType[region] == 'BED6': regionWorker(binFolder, region, chrom, chrToIndivRegions, limitSize, numBins, sideExtension, sideNumBins, readsForChrom, binLength)
		else: blockRegionWorker(binFolder, region, chrom, chrToIndivRegions, limitSize, numBins, sideExtension, sideNumBins, readsForChrom, binLength)
	
	logger.info('%s done: %s', chrom, ', '.join(regions.keys()))