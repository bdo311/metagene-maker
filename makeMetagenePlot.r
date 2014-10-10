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

# Extrapolate bins to exactly the number we need
# SNF: this function fails if there is only one bin currently, i.e. the line in test5 
extrap = function(x) {
  num = sum(!is.na(x))
  if (num==numBins) return(x)
  approx(1:num, x[1:num], n=numBins)$y
}

# extrap
no_col = count.fields(fn, sep = "\t") 
data = t(read.table(fn, sep='\t', fill=TRUE, col.names=1:max(no_col))) 

colnames(data) = data[nameCol,]
data = data[startCol:nrow(data),]
data = apply(data, c(1,2), as.numeric)
data.extrap = apply(data, 2, extrap) # coerce to exactly some number of bins

# Scale and collapse data, and print averages to files
data.avg = apply(data.extrap, 1, mean, na.rm=TRUE)
write.table(data.avg, paste('avgraw_', folderName, '_', regionName, '_', substr(fn, 1, nchar(fn)-4), ".txt", sep=''), sep='\t')

# Make plots of the scaled average data
pdf("plot.pdf")
plot(data.avg, type='l', lwd=2, ylab = "Average reads per nt", main="All regions", col="red")
dev.off()

save(data, data.extrap, data.avg, file="data.RData")
