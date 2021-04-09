#Code for RNAseq analysis of mouse placental tissue for the "Developmental Circadian Disruption Alters Placental Signaling in Mice" paper
#step 0: load libraries and data
setwd("your_path_here")
library(Hmisc)
library(gmodels)
library(DESeq2)
library(sva) 
library(ggplot2)
library(ggrepel)
library(dplyr)
#load("DESeq_mouse.RData")
#load gene count data
cts<-read.table("genecounts_placenta.txt", check.names = FALSE, sep="\t", header=TRUE) 
cts[1:2,]
head(cts,2)
rownames(cts) <- cts[,c(1)]
cts[1] <- NULL 
#load sample information
PHENO <- read.csv("Pheno_data_rnaseq.csv")
rownames(PHENO) <- PHENO$Sample_ID
save(PHENO, cts, file = "DESeq_mouse.RData")
#
load("DESeq_mouse.RData")
#
#step 1: format data for DESeq2 analysis
names(PHENO)
PHENO1 <- PHENO[ which(PHENO$Sample_ID!=1012), ] # drop the outlier sample, sample #12
PHENO1$Treatment <- as.factor(PHENO1$Treatment) #make sure factor variable
PHENO1$Sex <- as.factor(PHENO1$Sex)
# extract phenotypic variables and drop gene count data for the outlier sample
coldata <- PHENO1[,c("Sample_ID", "Treatment", "Sex", "Dam_ID", "Fetus_ID", "RIN", "Nucleic_conc", "s260_280", "s260_230", "Placental_wt", "Dam_wt", "sex_ratio_FtoM")]
coldata 
all(rownames(coldata) %in% colnames(cts)) #TRUE
all(rownames(coldata) == colnames(cts)) #FALSE
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts)) #TRUE
coldata
############################################################
#step 2: run the differential expression analysis
#specify model design
dds <-DESeqDataSetFromMatrix(countData = cts,
                             colData = coldata,
                             design = ~ Sex + Treatment) #put exposure variable of interest to far right
dds
# can test for interaction by changing design = ~ Sex + Treatment + Sex:Treatment
#QC filtering of RNA-seq data
isexpr <- rowSums(fpm(dds)>1) >= ncol(dds)*.10
dds <- dds[isexpr,] #filters # total transcripts from 53,801 to 14,743
#Evaluate results
dds$Sex <- relevel(dds$Sex, ref = "Male") #making "male" be the reference group
dds$Treatment <- relevel(dds$Treatment, ref = "CL") # making "CL" be reference group
#run the actual DE analysis with specified model design
dds <-DESeq(dds) #this is using Wald method so comparisons can be made
#
resultsNames(dds) #[1] "Intercept"          "Sex_Female_vs_Male" "Treatment_CD_vs_CL"
###############################
# step 3: run the model controlling for a surrogate variable
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~Sex + Treatment, colData(dds)) #full model with main exposure
mod0 <- model.matrix(~Sex, colData(dds)) #reduced model 
#use sva package to determine surrogate variables
n.sv=num.sv(dat=dat, mod=mod, method="be", B=200) #200 iterations to get more stable estimate, may take a few minutes
sva.seq = svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=n.sv, method="irw") #gives number of surrogate variables to be 2; to avoid overadjustment, we'll only use the first SV
#adding the first surrogate variable to the dataset for analysis
ddssva <- dds #create a new dataset with SVs included
ddssva$SV1 <- sva.seq$sv[,1] #first surrogate variable added to coldata
ddssva
#model now with SV
design(ddssva) <- ~ SV1 + Sex + Treatment
ddssva <- DESeq(ddssva) #run DESeq2
#
resultsNames(ddssva) #[1] "Intercept" "SV1" "Sex_Female_vs_Male" "Treatment_CD_vs_CL"
#######################################
#step 4: get the sex-specific DE gene results
resultsNames(ddssva)
res1 <- results(ddssva, name="Sex_Female_vs_Male") #Male is ref
#
resOrdered1 <- res1[order(res1$pvalue),] #order results in the table by smallest p-value
resOrdered1$padj
resOrdered1
summary(res1) #summarize basic tallies
sum(res1$padj < 0.05, na.rm=TRUE) #how many adj p-values <0.05?
#[1] 77
resOrdered1$padj_bon <- p.adjust(resOrdered1$pvalue, method = "bonferroni")#provides bonferroni-adjusted p-values
write.csv(as.data.frame(resOrdered1), file="Sex_Female_vs_Male.csv") #write results to file
#######################################
#step 5: get the light treatment-specific DE gene results
res <- results(ddssva, name="Treatment_CD_vs_CL") #CL is ref
res
#
resOrdered <- res[order(res$pvalue),] #order results in the table by smallest p-value
resOrdered$padj
resOrdered
summary(res) #summarize basic tallies
sum(res$padj < 0.05, na.rm=TRUE) #how many adj p-values <0.05?
#[1] 9
resOrdered$padj_bon <- p.adjust(resOrdered$pvalue, method = "bonferroni")
write.csv(as.data.frame(resOrdered), file="Treatment_CD_vs_CL.csv")
#
#######################################
#step 6: plot results
#PCA plot with treatment
rld <- rlog(ddssva, blind=FALSE)
z <- plotPCA(rld, intgroup="Treatment") #can change variable to color-code dots by
z

#volcano plot of sex-specific hits
res1 <- results(ddssva, name="Sex_Female_vs_Male") 
res3 <- as.data.frame(res1)
with(res3, plot(log2FoldChange, -log10(pvalue), pch=20,cex = 1, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col="#999999", xlim=c(-12, 5), ylim=c(0, 110), main="Male vs Female"))
with(subset(res3 , padj<.05), points(log2FoldChange, -log10(pvalue), pch=20,col="red")) #BH sig hits in red
with(subset(res3 , padj<0.002), points(log2FoldChange, -log10(pvalue), pch=20,col="black", cex=1.2)) #bonferroni sig hits in black

#volcano plot of treatment-specific hits
res <- results(ddssva, name="Treatment_CD_vs_CL") 
res2 <- as.data.frame(res)
with(res2, plot(log2FoldChange, -log10(pvalue), pch=20,cex = 1, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col="#999999", xlim=c(-3, 3), ylim=c(0, 10), main="CL vs CD"))
with(subset(res2 , padj<.05), points(log2FoldChange, -log10(pvalue), pch=20,col="red")) #BH sig hits in red
with(subset(res2 , padj<0.01), points(log2FoldChange, -log10(pvalue), pch=20,col="black", cex=1.2)) #bonferroni sig hits in black

#counts plot
topGene <- rownames(res)[which.min(res$padj)]
topGene = "ENSMUSG00000000753"
plotCounts(ddssva, gene = topGene, intgroup=c("Treatment"), main="Serpinf1",  xlab="Treatment",  cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex.names=1.5)

topGene = "ENSMUSG00000040489"
plotCounts(ddssva, gene = topGene, intgroup=c("Treatment"), main="Sox30", xlab="Treatment", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex.names=1.5 )
