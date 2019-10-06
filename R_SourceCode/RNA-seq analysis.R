if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("edgeR")
BiocManager::install("Mus.musculus")
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

#PART 4.1

#extracts the data from the article for analysis 
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

#Combining files into a large char vector with a column for each sample
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

#makes a HUGE matrix with 27,179 rows and 9 columns with all of the gene IDs and samples (respectively)
x <- readDGE(files, columns=c(1,3))
class(x)

#PART 4.2
#Adding experimental variables. In this case, adding cell type variables
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
#Adding the cell type variables
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

#PART 4.3
#Makes a 2nd data frame for gene-level info. Mus.musculus is for house mice, homo.sapiens for humans
#Adds gene symbols, names, chromosome names, locations etc
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)
#Removes duplicated genes, keeping only the first occurence of each gene
genes <- genes[!duplicated(genes$ENTREZID),]
x$genes <- genes
x

#PART 5.1
#Transforming the measurements into CPM and log-CPM
#CPM of 1 is lowest # counts, while a CPM of 76 is highest # counts
#Log-CPM used for plotting (complex math)/linear modeling in Limma
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#Getting the mean
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

#PART 5.2
#Removes low expressed genes. Includes genes that are not expressed or those that are very low expressed
#Biologically, this helps as low expressed genes are mostly irrelevant to the bodily processes
table(rowSums(x$counts==0)==9)
#Filters by most significant
#Uses CPM values, keeps genes with CPM of 0.2 or more expressed in 3 samples
#Lower CPM cutoff would have been used with a larger library, higher CPM cutoff with a smaller library
#Reduced to 60% of the original list of 16,624 genes
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#Plots the difference between the raw and unfiltered data
plot.new()
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")




#PART 5.3
#Values are normalised to ensure no external factors are measured 
#First batches often have higher overall expression than later tests
#Creates the normalization factors
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#Example of exaggerated normalization changes
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

#Plots the normalised and unormalized data
#UNORMALIZED
plot.new()
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
#NORMALIZED
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")


#PART 5.4
#Make a multi-dimensional scaling (MDS) plot to show diferences between samples
plot.new()
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

#PART 6.1
#Compares gene expression between the 3 cell populations profiled
#Assume data has normal distribution
#Design matrix created with cell population and sequencing lane (batch) info
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design
#Then make contrasts between the results using limma makeContrasts function (via linear modeling)
contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

#PART 6.2
#Variance independent of the mean, assume normal distribution
#Limma conerts raw counts to log-CPM values using normalisation factors from the data frame itself
#Creates a voom-plot, typically with decreasing trend between means and variances
#Voom-plot shows level of filtering. Drop in variance levels signifies poor filtering of low genes. If so, inicrease the dataset
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

#PART 6.4
#Summarizes the significant genes from the different cell populations
summary(decideTests(efit))
#Gets more significant genes through logFold Change measurement also
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
#makes a venn diagram
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))


#PART 6.5
#Most sigf genes are examined
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)
#PART 6.6
#Mean-difference plots (displaying logFC against average logCPM)
plot.new()
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
#Glimma makes an interactive plot, letting you search for specific genes
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)
#A heatmap is made to look at the expression of individual groups from a subset of genes. 
#It is made for the top 100 significant DE genes
plot.new()
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

