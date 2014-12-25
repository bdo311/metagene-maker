# combinePlusMinus.py
# 10/4/14
# combines + and - genes to get one plot for TSS, TES, and genebody

import csv, os, glob, collections
from multiprocessing import Pool

csv.register_dialect("textdialect", delimiter='\t')

# write the mapping to a file
def writeFile(name, mapping, direc):
	ofile = open(direc + name + '.txt', 'w')
	writer = csv.writer(ofile, 'textdialect')
	for treatment in mapping:
		outputRow = [treatment]
		outputRow.extend(mapping[treatment])
		writer.writerow(outputRow)

	ofile.close()

def processFile(fileName, isMinus):
	os.chdir("/home/raflynn/7SK/GROseq/combined")
	print fileName
	avgFile = glob.glob(fileName + "*.txt")[0]
	print avgFile
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

def processRegions(sample):
	regions = ["tss", "geneBody", "tes"]
	dir = "/home/raflynn/7SK/GROseq/"

	
	for region in regions:
		# split allchr into plus and minus
		plusdir = dir + sample + "_plus/bins/" + region + '/'
		os.chdir(plusdir)
		os.system("rm -f allchr_*.txt")
		os.system("awk -F '\t' '{print >> \"allchr_\" $3 \".txt\"}' allchr.txt")
		minusdir = dir + sample + "_minus/bins/" + region + '/'
		os.chdir(minusdir)
		os.system("rm -f allchr_*.txt")
		os.system("awk -F '\t' '{print >> \"allchr_\" $3 \".txt\"}' allchr.txt")
		
		# combine plus plus, minus minus for sense
		file1 = plusdir + "/allchr_+.txt"
		file2 = minusdir + "/allchr_-.txt"

		opath = dir + "combined/"
		ofile = sample + "_" + region + "_sense.txt"
		os.system("cat " + file1 + " " + file2 + " > " + opath + ofile)
		
		os.chdir(opath)
		cmd = "Rscript /home/raflynn/Scripts/metagene_maker/makeMetagenePlot.r " + sample + " " + region + " 6 0 100 " + ofile
		print cmd
		os.system(cmd)
		
		# combine plus minus, minus plus for antisense
		file1 = minusdir + "/allchr_+.txt"
		file2 = plusdir + "/allchr_-.txt"
		
		opath = dir + "combined/"
		ofile = sample + "_" + region + "_antisense.txt"
		os.system("cat " + file1 + " " + file2 + " > " + opath + ofile)
		
		os.chdir(opath)
		cmd = "Rscript /home/raflynn/Scripts/metagene_maker/makeMetagenePlot.r " + sample + " " + region + " 6 0 100 " + ofile
		print cmd
		os.system(cmd)
	
def main():
	samples = ["GRO_" + x + 'comb' for x in ["12C", "125", "123", "6C", "65", "63"]]
	#regions = ["cleavage500"]
	dir = "/home/raflynn/7SK/GROseq/"
	
	p = Pool(6)
	p.map(processRegions, samples)
	
	regions = ["tss", "geneBody", "tes"]
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	name = "groseq_asoComb_combined"
	for region in regions:
		for direc in ["_sense", "_antisense"]:
			os.chdir(dir)
			for s in samples:
				fn = "avgraw*_" + s + "_" + region + direc			
				folder = s + "_" + region + direc			
			
				if 'antisense' not in folder:
					regionToFolderAvgs[region][folder] = processFile(fn, isMinus=False)
				else: 
					regionToFolderAvgs[region][folder] = processFile(fn, isMinus=True)
		writeFile(name + '_' + region, regionToFolderAvgs[region], dir + '/averages/')
	
if __name__ == "__main__":
	main()