#!/usr/bin/env python

# Stephen N. Floor 
# 7 October 2014 
# floor@berkeley.edu 

# TODO: 

# convert UCSC gene names to refseq?
# add a region type which is mrna that contains the whole spliced transcript, and preserve cdsStart and cdsEnd for these.  - this is the same as codingExons, just need to preserve start/stop 

import sys, os, csv
from UCSCKnownGene import * 
csv.register_dialect("textdialect", delimiter='\t') #BTD added for later parsing

def main():
	print " ---------------------------------"
	print "| Extract Regions from knownGenes |"
	print "|  snf   7 October 2014           |"
	print " ---------------------------------\n\n"

	if len(sys.argv) != 3:
		sys.exit("ERROR: Please provide the knownGenes file and output basename as inputs")


	# output filenames: 
	utr5FName = sys.argv[2] + "_5utr.bed"
	utr5StartFName = sys.argv[2] + "_5utr_start.bed"
	cdsFName = sys.argv[2] + "_cds.bed"
	utr3FName = sys.argv[2] + "_3utr.bed"
	exonFName = sys.argv[2] + "_exons.bed"
	intronFName = sys.argv[2] + "_introns.bed"
	codingExonFName = sys.argv[2] + "_codingexons.bed"
	codingIntronFName = sys.argv[2] + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA 
	noncodingExonFName = sys.argv[2] + "_noncodingexons.bed" 
	noncodingIntronFName = sys.argv[2] + "_noncodingintrons.bed" 

	#keep track of where we are 
	genesRead = 0

	# parameters that should be passed via the cmd line
	useBlocks = True 


	# #terminate if output files exist

	# if os.path.exists(utr5FName) or os.path.exists(utr5StartFName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) \
			# or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
		# sys.exit("ERROR: output basename %s files already exist" % sys.argv[2]) 

	# #process the file

	with open(sys.argv[1]) as knownGenesFile, open(utr5FName, "w") as utr5File, open(utr5StartFName, "w") as utr5StartFile, open(cdsFName, "w") as cdsFile, \
			open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
			open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
			open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile:

		for line in knownGenesFile:
			# all of the knowngenes parsing and metadata construction is done inside UCSCKnownGene.py, especially the createGene method 
			gene = createGene(line)
			genesRead += 1

			if (useBlocks): # output all region primitives on the same line by specifying nBlocks and lists inside the BED output 
				if(gene.coding): 
					#blockBedFormat is one line by definition 
					if (gene.utr5Len > 0): utr5File.write(gene.blockBedFormat(region="5utr") + "\n")
					if (gene.utr5startLen > 0): utr5StartFile.write(gene.blockBedFormat(region="5utr_start") + "\n")
					if (gene.cdsLen > 0): cdsFile.write(gene.blockBedFormat(region="cds") + "\n")
					if (gene.utr3Len > 0): utr3File.write(gene.blockBedFormat(region="3utr") + "\n")
					
					if (gene.exonsLen > 0): 
						exonFile.write(gene.blockBedFormat(region="exons") + "\n")
						codingExonFile.write(gene.blockBedFormat(region="exons") + "\n")
					
					if (gene.intronsLen > 0):
						intronFile.write(gene.blockBedFormat(region="introns") + "\n")
						codingIntronFile.write(gene.blockBedFormat(region="introns") + "\n")
				
				else: # noncoding transcripts just have exons and introns 
					if (gene.exonsLen > 0):
						exonFile.write(gene.blockBedFormat(region="exons") + "\n")
						noncodingExonFile.write(gene.blockBedFormat(region="exons") + "\n")
					
					if (gene.intronsLen > 0):
						intronFile.write(gene.blockBedFormat(region="introns") + "\n")
						noncodingIntronFile.write(gene.blockBedFormat(region="introns") + "\n")

			else: # output one line per region primitive instead of combining regions via blocks 
				if(gene.coding): 
					for entry in gene.bedFormat(region="5utr"):
						utr5File.write(entry + "\n")
					for entry in gene.bedFormat(region="5utr_start"):
						utr5StartFile.write(entry + "\n")
					for entry in gene.bedFormat(region="cds"):
						cdsFile.write(entry + "\n")
					for entry in gene.bedFormat(region="3utr"):
						utr3File.write(entry + "\n")

					for entry in gene.bedFormat(region="exons"):
						exonFile.write(entry + "\n")
						codingExonFile.write(entry + "\n")
				
					for entry in gene.bedFormat(region="introns"):
						intronFile.write(entry + "\n")
						codingIntronFile.write(entry + "\n")
				
				else: # noncoding transcripts just have exons and introns 
					for entry in gene.bedFormat(region="exons"):
						exonFile.write(entry + "\n")
						noncodingExonFile.write(entry + "\n")
				
					for entry in gene.bedFormat(region="introns"):
						intronFile.write(entry + "\n")
						noncodingIntronFile.write(entry + "\n")

			if (not genesRead % 2500):
				print "Processed %d entries..." %  genesRead

	print "Processed %d entries." %  genesRead

	# BTD Edit: making unique regions and linking to gene name
	# --------------------------------------------------------
	# utr5FName = sys.argv[2] + "_5utr.bed"
	# utr5StartFName = sys.argv[2] + "_5utr_start.bed"
	# cdsFName = sys.argv[2] + "_cds.bed"
	# utr3FName = sys.argv[2] + "_3utr.bed"
	# exonFName = sys.argv[2] + "_exons.bed"
	# intronFName = sys.argv[2] + "_introns.bed"
	# codingExonFName = sys.argv[2] + "_codingexons.bed"
	# codingIntronFName = sys.argv[2] + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA 
	# noncodingExonFName = sys.argv[2] + "_noncodingexons.bed" 
	# noncodingIntronFName = sys.argv[2] + "_noncodingintrons.bed" 

	# 1. Get gene ID (NM_123, ENSG123) --> gene name (Abcd1)
	print "Getting gene ID"
	idToName = {}
	with open(sys.argv[1], 'r') as knownGeneFile:
		reader = csv.reader(knownGeneFile, 'textdialect')
		for row in reader:
			idToName[row[0]] = row[-1]
			
	# 2. Get unique 5'UTR, 5'Start UTR, and 3'UTR
	print "Getting unique UTRs"
	def getUniqUTR(uniqFN, utrFN):
		with open(uniqFN, 'w') as uniq, open(utrFN, 'r') as utr:
			already = set()
			reader = csv.reader(utr, 'textdialect')
			writer = csv.writer(uniq, 'textdialect')
			for row in reader:
				if tuple(row[6:]) in already: continue #repeat
				geneIDInfo = row[3]
				id = geneIDInfo.split('__')[0]
				geneName = idToName[id]
				row[3] = geneIDInfo + '__' + geneName
				already.add(tuple(row[6:]))
				writer.writerow(row)
				
	uniq5UTR = sys.argv[2] + "_uniq_5utr.bed"
	getUniqUTR(uniq5UTR, utr5FName)

	uniq3UTR = sys.argv[2] + '_uniq_3utr.bed'
	getUniqUTR(uniq3UTR, utr3FName)

	uniq5SUTR = sys.argv[2] + '_uniq_5utr_start.bed'
	getUniqUTR(uniq5SUTR, utr5StartFName)
		
	# 3. Get unique exons + num. Do it 3x for all, coding, and noncoding
	print "Getting unique exons"
	def getUniqExons(uniqFN, exonFN):
		with open(uniqFN, 'w') as uniq, open(exonFN, 'r') as exons:
			already = set()
			reader = csv.reader(exons, 'textdialect')
			writer = csv.writer(uniq, 'textdialect')
			for row in reader:
				# gene ID info
				geneIDInfo = row[3]
				id = geneIDInfo.split('__')[0]
				geneName = idToName[id]
				geneIDInfo = geneIDInfo + '__' + geneName
				
				# chrom, start, stop, strand
				chrom = row[0]
				start, end = int(row[1]), int(row[2])
				strand = row[5]

				# calculate exon starts and lengths
				exonLengths = row[10].split(',')
				if exonLengths[-1] == '': exonLengths = exonLengths[:-1]
				exonLengths = [int(x) for x in exonLengths]
				exonStarts = row[11].split(',')
				if exonStarts[-1] == '': exonStarts = exonStarts[:-1]
				exonStarts = [int(x) for x in exonStarts]
				
				# calculate exons
				exons = []
				for i in range(len(exonStarts)):
					absStart = start + exonStarts[i]
					exons.append([absStart, absStart + exonLengths[i]])
				if strand == '-': exons = exons[::-1] #flip exon order
				
				# making BED6
				for i in range(len(exons)):
					exonNum = i + 1
					exon = exons[i]
					outputRow = [chrom, exon[0], exon[1]]
					
					# unique
					if tuple(outputRow) in already: continue
					already.add(tuple(outputRow))			
					outputRow.extend([geneIDInfo, exonNum, strand])
					writer.writerow(outputRow)
				
	uniqExons = sys.argv[2] + '_uniq_exons.bed'
	getUniqExons(uniqExons, exonFName)

	uniqExons = sys.argv[2] + '_uniq_codingexons.bed'
	getUniqExons(uniqExons, codingExonFName)

	uniqExons = sys.argv[2] + '_uniq_noncodingexons.bed'
	getUniqExons(uniqExons, noncodingExonFName)			

	# 4. Get unique introns + num. unique 5'SS, 3'SS. 
	# 5'SS is first base of intron, 3'SS is last base of intron
	print "Getting unique introns and 5' and 3' SS"
	def getUniqIntronsAndSS(uniqIntronFN, uniq5SSFN, uniq3SSFN, intronFN):
		with open(uniqIntronFN, 'w') as uniqIntron, open(uniq5SSFN, 'w') as uniq5, \
			open(uniq3SSFN, 'w') as uniq3, open(intronFN, 'r') as introns:
			alreadyIntron = set()
			already5 = set()
			already3 = set()
			
			reader = csv.reader(introns, 'textdialect')
			intronWriter = csv.writer(uniqIntron, 'textdialect')
			fiveWriter = csv.writer(uniq5, 'textdialect')
			threeWriter = csv.writer(uniq3, 'textdialect')
			
			for row in reader:
				# gene ID info
				geneIDInfo = row[3]
				id = geneIDInfo.split('__')[0]
				geneName = idToName[id]
				geneIDInfo = geneIDInfo + '__' + geneName
				
				# chrom, start, stop, strand
				chrom = row[0]
				start, end = int(row[1]), int(row[2])
				strand = row[5]

				# calculate intron starts and lengths
				intronLengths = row[10].split(',')
				if intronLengths[-1] == '': intronLengths = intronLengths[:-1]
				intronLengths = [int(x) for x in intronLengths]
				intronStarts = row[11].split(',')
				if intronStarts[-1] == '': intronStarts = intronStarts[:-1]
				intronStarts = [int(x) for x in intronStarts]
				
				# calculate introns
				introns = []
				for i in range(len(intronStarts)):
					absStart = start + intronStarts[i]
					introns.append([absStart, absStart + intronLengths[i]])
				if strand == '-': introns = introns[::-1] #flip intron order
				
				# making BED6
				for i in range(len(introns)):
					intronNum = i + 1
					intron = introns[i]
					outputRow = [chrom, intron[0], intron[1]]
					
					# unique introns
					if tuple(outputRow) in alreadyIntron: continue
					alreadyIntron.add(tuple(outputRow))
					outputRow.extend([geneIDInfo, intronNum, strand])
					intronWriter.writerow(outputRow)
					
					# unique splice sites
					if strand == '+':
						fiveSS = [chrom, intron[0], intron[0] + 1]
						threeSS = [chrom, intron[1] - 1, intron[1]]
					else:
						threeSS = [chrom, intron[0], intron[0] + 1]
						fiveSS = [chrom, intron[1] - 1, intron[1]]
					if tuple(fiveSS) not in already5:
						already5.add(tuple(fiveSS))
						fiveSS.extend([geneIDInfo, intronNum, strand])
						fiveWriter.writerow(fiveSS)
					if tuple(threeSS) not in already3:
						already3.add(tuple(threeSS))
						threeSS.extend([geneIDInfo, intronNum, strand])
						threeWriter.writerow(threeSS)

	uniqIntrons = sys.argv[2] + '_uniq_introns.bed'
	uniq5 = sys.argv[2] + '_uniq_5ss.bed'
	uniq3 = sys.argv[2] + '_uniq_3ss.bed'
	getUniqIntronsAndSS(uniqIntrons, uniq5, uniq3, intronFName)

	uniqIntrons = sys.argv[2] + '_uniq_codingintrons.bed'
	uniq5 = sys.argv[2] + '_uniq_coding5ss.bed'
	uniq3 = sys.argv[2] + '_uniq_coding3ss.bed'
	getUniqIntronsAndSS(uniqIntrons, uniq5, uniq3, codingIntronFName)

	uniqIntrons = sys.argv[2] + '_uniq_noncodingintrons.bed'
	uniq5 = sys.argv[2] + '_uniq_noncoding5ss.bed'
	uniq3 = sys.argv[2] + '_uniq_noncoding3ss.bed'
	getUniqIntronsAndSS(uniqIntrons, uniq5, uniq3, noncodingIntronFName)

	# 5. unique TSS/TES
	print "Getting unique TSS and TES"
	def getUniqTSSAndTES(tssFN, tesFN, cdsFN):
		with open(tssFN, 'w') as uniqTSS, open(tesFN, 'w') as uniqTES, open(cdsFN, 'r') as cds:
			alreadyTSS = set()
			alreadyTES = set()
			reader = csv.reader(cds, 'textdialect')
			tssWriter = csv.writer(uniqTSS, 'textdialect')
			tesWriter = csv.writer(uniqTES, 'textdialect')
			for row in reader:
				geneIDInfo = row[3]
				id = geneIDInfo.split('__')[0]
				geneName = idToName[id]
				geneIDInfo = geneIDInfo + '__' + geneName
				
				# chrom, start, stop, strand
				chrom = row[0]
				strand = row[5]
				start, end = int(row[1]), int(row[2])
				
				if strand == '+':
					startRow = [chrom, start, start + 1]
					endRow = [chrom, end - 1, end]
				else:
					startRow = [chrom, end - 1, end]
					endRow = [chrom, start, start + 1]
				if tuple(startRow) not in alreadyTSS:
					alreadyTSS.add(tuple(startRow))
					startRow.extend([geneIDInfo, 0, strand])
					tssWriter.writerow(startRow)
				if tuple(endRow) not in alreadyTSS:
					alreadyTES.add(tuple(endRow))
					endRow.extend([geneIDInfo, 0, strand])
					tesWriter.writerow(endRow)			
				
	uniqTSS = sys.argv[2] + '_uniq_tss.bed'
	uniqTES = sys.argv[2] + '_uniq_tes.bed'
	getUniqTSSAndTES(uniqTSS, uniqTES, cdsFName)
	
if __name__ == '__main__':
	main()
	