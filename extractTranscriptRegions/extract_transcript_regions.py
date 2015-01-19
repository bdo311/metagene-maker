#!/usr/bin/env python

# Stephen N. Floor
# 7 October 2014
# floor@berkeley.edu

# TODO:

# convert UCSC gene names to refseq?
# add a region type which is mrna that contains the whole spliced transcript, and preserve cdsStart and cdsEnd for these.  - this is the same as codingExons, just need to preserve start/stop
# remove pandas from GTF.py 

import sys, os, argparse, csv, glob
csv.register_dialect("textdialect", delimiter='\t') #BTD added for later parsing
import GTF
from collections import defaultdict 
from Transcript import *

print " ----------------------------------"
print "| Extract Regions from annotations |"
print "|  snf        Fall 2014            |"
print " ----------------------------------\n\n"


# ------ ARGUMENT PARSING ----------

parser = argparse.ArgumentParser(description="Create transcript regions (5' UTR/CDS/3'UTR etc) from knownGenes or a GTF") 

parser.add_argument("-i", "--input", help="input filename", required=True)
parser.add_argument("-o", "--output", help="output basename", required=True) 
parser.add_argument("--ucsc", help="Read from a UCSC knownGenes formatted file (BED)", action="store_true")
parser.add_argument("--gtf", help="Read from a GTF (only tested with Ensembl GTFs)", action="store_true")

args = parser.parse_args() 

if ( (not (args.ucsc or args.gtf)) or (args.ucsc and args.gtf)):
    sys.exit("FATAL: must set one but not both of --ucsc and --gtf") 

# output filenames:
utr5FName = args.output + "_5utr.bed"
utr5StartFName = args.output + "_5utr_start.bed"
cdsFName = args.output + "_cds.bed"
utr3FName = args.output + "_3utr.bed"
exonFName = args.output + "_exons.bed"
intronFName = args.output + "_introns.bed"
codingExonFName = args.output + "_codingexons.bed"
codingIntronFName = args.output + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA
noncodingExonFName = args.output + "_noncodingexons.bed"
noncodingIntronFName = args.output + "_noncodingintrons.bed"

#keep track of where we are

# parameters that should be passed via the cmd line
useBlocks = True


# #terminate if output files exist

# if os.path.exists(utr5FName) or os.path.exists(utr5StartFName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) \
        # or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
    # sys.exit("ERROR: output basename %s files already exist" % args.output)

# #process the file

def main():
	with open(utr5FName, "w") as utr5File, open(utr5StartFName, "w") as utr5StartFile, open(cdsFName, "w") as cdsFile, \
			open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
			open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
			open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile:

		def writeOutput(gene):
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


		if (args.ucsc): 
			with open(args.input, "r") as genesFile: 
				genesRead = 0

				for line in genesFile:
					# all of the knowngenes parsing and metadata construction is done inside UCSCKnownGene.py, especially the createGene method

					gene = createUCSCTranscript(line) 
					genesRead += 1

					writeOutput(gene)

					if (not genesRead % 2500):
						print "Processed %d entries..." %  genesRead

					
		elif (args.gtf): 
				
				# first parse the entire file into a dictionary of lists

			txDict = defaultdict(list) 

			print "Building GTF dictionary..." 

			# the issue here is that lines for various transcripts may be interleaved, so can either create lots of SNFGene objects, or a giant dict. opted for giant dict. 
			for line in GTF.lines(args.input): 

				txDict[line["transcript_id"]].append(line)
				genesRead += 1

				if (not genesRead % 100000):
					print "Processed %d lines..." %  genesRead

			print "Dictionary built." 

			# now create a SNFGene object for each transcript and output it 
			genesRead = 0
			for key in txDict: 

				#print key

				tx = createGTFTranscript(txDict[key])

				#print tx 
				writeOutput(tx)
				genesRead += 1
				
				if (not genesRead % 2500):
					print "Processed %d entries..." %  genesRead


	print "Processed %d entries." %  genesRead

	# BTD Edit: making unique regions and linking to gene name
	# --------------------------------------------------------
	# utr5FName = args.output  + "_5utr.bed"
	# utr5StartFName = args.output  + "_5utr_start.bed"
	# cdsFName = args.output  + "_cds.bed"
	# utr3FName = args.output  + "_3utr.bed"
	# exonFName = args.output  + "_exons.bed"
	# intronFName = args.output  + "_introns.bed"
	# codingExonFName = args.output  + "_codingexons.bed"
	# codingIntronFName = args.output  + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA 
	# noncodingExonFName = args.output  + "_noncodingexons.bed" 
	# noncodingIntronFName = args.output  + "_noncodingintrons.bed" 

	# 1. Get gene ID (NM_123, ENSG123) --> gene name (Abcd1)
	print "Getting gene ID"
	idToName = {}
	if args.ucsc:
		with open(args.input, 'r') as knownGeneFile:
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
				try: geneName = idToName[id]
				except: geneName = id
				if geneName != id: row[3] = id + '__' + geneName
				else: row[3] = id
				already.add(tuple(row[6:]))
				writer.writerow(row)
				
	uniq5UTR = args.output  + "_uniq_5utr.bed"
	getUniqUTR(uniq5UTR, utr5FName)

	uniq3UTR = args.output  + '_uniq_3utr.bed'
	getUniqUTR(uniq3UTR, utr3FName)

	uniq5SUTR = args.output  + '_uniq_5utr_start.bed'
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
				try: geneName = idToName[id]
				except: geneName = id
				if geneName != id: geneIDInfo = id + '__' + geneName
				else: geneIDInfo = id
				
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
					exonNumInfo = str(exonNum) + 'of' + str(len(exons))
					exon = exons[i]
					outputRow = [chrom, exon[0], exon[1]]
					
					# unique
					if tuple(outputRow) in already: continue
					already.add(tuple(outputRow))            
					outputRow.extend([geneIDInfo + '__exon__' + exonNumInfo, 0, strand])
					writer.writerow(outputRow)
				
	uniqExons = args.output  + '_uniq_exons.bed'
	getUniqExons(uniqExons, exonFName)

	uniqExons = args.output  + '_uniq_codingexons.bed'
	getUniqExons(uniqExons, codingExonFName)

	uniqExons = args.output  + '_uniq_noncodingexons.bed'
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
				try: geneName = idToName[id]
				except: geneName = id
				if geneName != id: geneIDInfo = id + '__' + geneName
				else: geneIDInfo = id
				
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
					intronNumInfo = str(intronNum) + 'of' + str(len(introns))
					intron = introns[i]
					outputRow = [chrom, intron[0], intron[1]]
					
					# unique introns
					if tuple(outputRow) in alreadyIntron: continue
					alreadyIntron.add(tuple(outputRow))
					outputRow.extend([geneIDInfo+ '__intron__' + intronNumInfo, 0, strand])
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
						fiveSS.extend([geneIDInfo + '__5ss__' + intronNumInfo, 0, strand])
						fiveWriter.writerow(fiveSS)
					if tuple(threeSS) not in already3:
						already3.add(tuple(threeSS))
						threeSS.extend([geneIDInfo+ '__3ss__' + intronNumInfo, 0, strand])
						threeWriter.writerow(threeSS)

	uniqIntrons = args.output  + '_uniq_introns.bed'
	uniq5 = args.output  + '_uniq_5ss.bed'
	uniq3 = args.output  + '_uniq_3ss.bed'
	getUniqIntronsAndSS(uniqIntrons, uniq5, uniq3, intronFName)

	uniqIntrons = args.output  + '_uniq_codingintrons.bed'
	uniq5 = args.output  + '_uniq_coding5ss.bed'
	uniq3 = args.output  + '_uniq_coding3ss.bed'
	getUniqIntronsAndSS(uniqIntrons, uniq5, uniq3, codingIntronFName)

	uniqIntrons = args.output  + '_uniq_noncodingintrons.bed'
	uniq5 = args.output  + '_uniq_noncoding5ss.bed'
	uniq3 = args.output  + '_uniq_noncoding3ss.bed'
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
				try: geneName = idToName[id]
				except: geneName = id
				if geneName != id: geneIDInfo = id + '__' + geneName
				else: geneIDInfo = id
				
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
				
	uniqTSS = args.output  + '_uniq_tss.bed'
	uniqTES = args.output  + '_uniq_tes.bed'
	getUniqTSSAndTES(uniqTSS, uniqTES, cdsFName)


	# sort everything
	print "Sorting BED files"
	for fn in glob.glob("*.bed"):
		os.system("sort -k1,1 -k2,2n %s -o %s"%(fn, fn))
		
if __name__=='__main__':
	main()