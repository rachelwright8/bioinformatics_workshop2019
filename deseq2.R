# Install the bioconductor packages (you only need to do this once!)
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2","arrayQualityMetrics"))

library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers

# Install these from CRAN, if you haven't already
library("pheatmap") # for pretty heatmaps
library("tidyverse") # for wrangling, plotting
library("VennDiagram") # for making Venn diagram
library("vegan") # for distance matrix calculations
library("ape") # for PCoA

setwd("/Volumes/wrightworkshop/rachel/") # change to your working directory

# Load data ---------------------------------------------------------------
countdata <- read.table("../shared_files/countsTable.csv",header=TRUE,row.names = 1) 
head(countdata) 
length(countdata[,1])
# 62559 isogroups ("genes") mapped

# simplify sample names
names(countdata)
names(countdata) <- gsub("MC","",names(countdata)) # replace the "MC" with nothing
names(countdata) <- gsub(".fq.trim.sam.counts","",names(countdata)) # replace the ".fq.trim.sam.counts" with nothing
names(countdata)

# make conditions table in excel and load it here
# or make it here if it's a simple design
conds <- read.csv("../shared_files/conds.csv")
head(conds)

# VERY IMPORTANT!
# The sample names in the "conds" table must match the sample names in the counts matrix (column names)
# The order must be the same. Be sure to check it!

# First, check that the names in the countdata file exist as sample names in conds
table(names(countdata) %in% conds$sam)
# Yes! All of the names in count data are also in the conds table

# Next, check that the order matches
table(names(countdata) == conds$sam)
# No! The order is totally off. Let's fix it.

# Reorder the samples conds table to match the order in the counts matrix
conds <- conds[match(names(countdata), conds$sam),]
head(conds)

# check that it worked
table(names(countdata) == conds$sam)
# yes!

# now that everything matches, lets make our lives easier by giving our samples more informative names
head(conds)
conds$name <- paste(conds$sam,conds$bank,conds$pheno,sep="_")
head(conds)

names(countdata) <- conds$name
head(countdata)

# ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts <- colSums(countdata)
min(totalCounts) # 126030
mean(totalCounts) # 568377.5
max(totalCounts)  # 956429

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = conds,
  design = ~ pheno + bank)

save(conds, countdata, dds, file="ddsCoral_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(countdata,1,mean)
table(means>3)
# FALSE  TRUE 
# 48052 14507  

means3 <- names(means[means>3])
head(means3)
length(means3)
#14507

coral.countFilt <- countdata[row.names(countdata) %in% means3,]
head(coral.countFilt)

totalCountsFilt <- colSums(coral.countFilt)
totalCountsFilt

min(totalCountsFilt) #122277
max(totalCountsFilt) #914705
mean(totalCountsFilt) #546563.9

# check sample order one more time... just in case!
table(names(coral.countFilt) == as.vector(conds$name))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = coral.countFilt,
  colData = conds,
  design = ~ pheno + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("pheno", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out2 <- c("11A_e_s", "17A_w_s", "18B_w_s", "1B_w_s", "3A_w_r", "7A_e_s")

# remove out3 and out2
dim(coral.countFilt)
countdata.out <- coral.countFilt %>% select(-one_of(c(out2)))
dim(countdata.out)

dim(conds)
head(conds)
conds.out <- conds %>% filter(!name %in% c(out2))
dim(conds.out)

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = countdata.out,
  colData = conds.out,
  design =  ~ pheno + bank)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)

#---Results
resPheno <- results(deds, independentFiltering = F, contrast=c("pheno","r","s"))
resPheno
# log2 fold change (MLE): pheno r vs s 
# Wald test p-value: pheno r vs s 
# here "resistant" is the numerator and "susceptible" is the denominator
# in other words, genes with a positive log2FoldChange are more highly expressed in resistant compared to susceptible genets

resBank <- results(deds, independentFiltering = F, contrast=c("bank", "e", "w"))
resBank
# log2 fold change (MLE): bank e vs w 
# Wald test p-value: bank e vs w 
# here "east" is the numerator and "west" is the denominator
# in other words, genes with a positive log2FoldChange are more highly expressed in east bank compared to west bank genets

# make a new VSD file (this is a useful format for making heatmaps later)
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# Write results for making heatmaps and other downstream analyses ---------------------------------------

###--------------Get pvals
head(resPheno)
valsPheno <- cbind(resPheno$pvalue, resPheno$padj)
head(valsPheno)
colnames(valsPheno)=c("pval.pheno", "padj.pheno")


head(resBank)
valsBank <- cbind(resBank$pvalue, resBank$padj)
head(valsBank)
colnames(valsBank)=c("pval.bank", "padj.bank")

#Make rlogdata and pvals table
vsdpvals <- as.data.frame(cbind(assay(vsd),valsPheno, valsBank))
head(vsdpvals)
dim(vsdpvals)
# 14507    52


# SAVE and LOAD (start here for downstream analyses) -------
save(conds.out, countdata.out, ddsFiltOut, deds,
     resPheno, resBank, vsd, vsdpvals, file="ddsCoralFiltOut.Rdata")

load("ddsCoralFiltOut.Rdata")

# Look at the results -------

# how many DE genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resPheno$padj < 0.05)
# FALSE  TRUE 
# 14486    21 

table(resBank$padj < 0.05)
# FALSE  TRUE 
# 14399   108 

# how many DE genes have log2 fold changes > 2 in either direction?
table(abs(resPheno$log2FoldChange)>1.5)
# FALSE  TRUE 
# 14429    78 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE 
# 14420    87 

# Explore with plots ------------------------------------------------------

#Sample distance heatmap
pheatmap(cor(assay(vsd)),border_color=NA, main="SampleHeatmap")

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds, main="Dispersion Plot")

#MA plot
plotMA(resPheno, ylim = c(-2, 2), main="MA Plot Pheno") 
plotMA(resBank, ylim = c(-2, 2), main="MA Plot Bank") 

# Venn diagram ---------
sig.pheno <-  row.names(vsdpvals[vsdpvals$padj.pheno<0.05 & !is.na(vsdpvals$padj.pheno),])
sig.bank <- row.names(vsdpvals[vsdpvals$padj.bank<0.05 & !is.na(vsdpvals$padj.bank),])

candidates <- list("Phenotype"=sig.pheno,"Bank"=sig.bank)

prettyvenn <- venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("blue", "forestgreen"),
  alpha = 0.5,
  label.col = c("black", "black", "black"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("black", "black"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)


# Make a gene expression heatmaps ------------
load("ddsCoralFiltOut.Rdata")

head(vsdpvals)
names(vsdpvals)
# subset for just the expression data
exp <- vsdpvals[c(1:48)]
head(exp)

#     Make p-value cut-offs
sig.pheno <-  row.names(vsdpvals[vsdpvals$padj.pheno<0.05 & !is.na(vsdpvals$padj.pheno),])
sig.bank <- row.names(vsdpvals[vsdpvals$padj.bank<0.01 & !is.na(vsdpvals$padj.bank),])

# subset expression data for significant DEGs
exp.pheno <- exp[row.names(exp) %in% sig.pheno,]
nrow(exp.pheno)

exp.bank <- exp[row.names(exp) %in% sig.bank,]
nrow(exp.bank)

# phenotype heatmap ---------
means <- apply(exp.pheno,1,mean) # means of rows
explc <- exp.pheno-means # subtracting them

# make colors 
heat.colors <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.1)(100)

annot_col <- data.frame(Phenotype = conds.out$pheno)
rownames(annot_col) <- names(explc)
annot_color <- list(Phenotype = c(r = "grey", s="black"))

# cluster plot
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1],
         clustering_distance_rows="correlation", cluster_cols=T)


# bank heatmap ---------
means <- apply(exp.bank,1,mean) # means of rows
explc <- exp.bank-means # subtracting them

# make colors 
heat.colors <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.6)(100)

annot_col <- data.frame(Bank = conds.out$bank)
rownames(annot_col) <- names(explc)
annot_color <- list(Bank = c(e = "grey", w="black"))


# cluster plot
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1],
         clustering_distance_rows="correlation", cluster_cols=T)

# PCoA ---------
load("ddsCoralFiltOut.Rdata")

# get variance stabilized expression data
head(vsdpvals)
names(vsdpvals)
# subset for just the expression data
exp <- vsdpvals[c(1:48)]
head(exp)

# make sure condition data match expression data
table(conds.out$name == names(exp))

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis(t(exp)~pheno+bank,
                    data=conds.out,method="manhattan")
adonisRes
#           Df  SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)   
# pheno      1   33144637 33144637  1.3358 0.02785    0.024 * 
# bank       1   40207141 40207141  1.6204 0.03379    0.002 **
# Residuals 45 1116586659 24813037         0.93836          
# Total     47 1189938438                  1.00000 

# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors

# plotting PCoA----
margin <- 1.5

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA for mid by site type
plot(scores[,xaxis], scores[,yaxis],type="n", 
     main = "Coral Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.axis=1.5,
     cex.main=1.5,
     cex.lab=1.5)
# plot "spiders" connecting the samples from each bank
ordispider(scores,conds.out$bank,label=T,cex=1.5)
# plot the resistant samples
points(scores[conds.out$pheno=="r",xaxis],
       scores[conds.out$pheno=="r",yaxis],
       col="grey", pch=15, cex=1.5) +
# plot the susceptible samples
points(scores[conds.out$pheno=="s",xaxis],
         scores[conds.out$pheno=="s",yaxis],
         col="black", pch=15, cex=1.5)


# legend of sites 
legend("topright", 
       c("Resistant", "Susceptible"),
       pch=c(15,15), 
       col=c("grey","black"), cex=1.5, bty = "n")

#insert p value 
legend("topleft",
       paste("Phenotype p = ",adonisRes$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=1.5, bty='n')  

legend("bottomleft", 
       paste("Bank p = ",adonisRes$aov.tab$`Pr(>F)`[2], sep=" "), 
       cex=1.5, bty='n')
