# File: makeMetagenePlot.r
# 8/29/14
# In a given folder, produces plot files and saves data into an RData

args = commandArgs(TRUE)
folderName = args[1]
regionName = args[2]
startCol = strtoi(args[3]) + 2
nameCol = strtoi(args[4]) + 1
numBins = strtoi(args[5])
fn = args[6]

extrap = function(x) {
  num = sum(!is.na(x))
  if (num==numBins) return(x)
  approx(1:num, x[1:num], n=numBins)$y
}

# extrap
data = read.table(fn, sep='\t', fill=TRUE, header=FALSE, colClasses=c(rep("character",6),rep("numeric",numBins+1)))

rownames(data) = make.names(data[,nameCol], unique=TRUE)
data = data[,startCol:ncol(data)]

if (numBins > 1) {
	data.extrap = t(apply(data, 1, extrap)) # coerce to exactly some number of bins

	# Scale and collapse data, and print averages to files
	data.avg = apply(data.extrap, 2, mean, na.rm=TRUE)
	write.table(data.avg, paste('avgraw_', folderName, '_', regionName, '_', substr(fn, 1, nchar(fn)-4), ".txt", sep=''), sep='\t', quote=FALSE)

	# Make plots of the scaled average data
	pdf("plot.pdf")
	plot(data.avg, type='l', lwd=2, ylab = "Average reads per nt", main="All regions", col="red")
	dev.off()
} else {
	data.avg = mean(data)
	write.table(data.avg, paste('avgraw_', folderName, '_', regionName, '_', substr(fn, 1, nchar(fn)-4), ".txt", sep=''), sep='\t', quote=FALSE)
}

save(data, data.extrap, data.avg, file="data.RData")
