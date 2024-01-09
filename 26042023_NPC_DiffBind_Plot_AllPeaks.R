library(DiffBind)

#Set working directory
setwd("~/Library/CloudStorage/OneDrive-Menntaský/Documents/ColdResponse-Salvör/CUT&RUN_analysis_in_R_NPCs")

#Reading in peaksets

NPC_H3K36_bg <- paste0("Bedgraph files/",list.files("Bedgraph files/", pattern = "*H3K36*"))
bamReads <- paste0("Bam files/NPC_H3K36_bam/",grep(list.files(path="Bam files/NPC_H3K36_bam/"), pattern='bai', invert=TRUE, value=TRUE))
SampleID <- substr(list.files("Bedgraph files/", pattern = "*H3K36*"), 1, 21)
Peaks <- paste0("NPC_SEACR_bed_files/Non_norm/", intersect(list.files("NPC_SEACR_bed_files/Non_norm/", pattern = "H3K36"), 
                                                           list.files("NPC_SEACR_bed_files/Non_norm/", pattern = "rep1|rep2")))

sample.sheet <- as.data.frame(cbind(SampleID,bamReads, Peaks))
sample.sheet$Condition <- c(32, 37, 32, 37)
sample.sheet$Replicate <- c(1,1,2,2)
sample.sheet$PeakCaller <- rep(c("bed"), times = 4)
sample.sheet$bamControl <- paste0("Bam files/NPC_IgG_bam/", grep(list.files(path = "Bam files/NPC_IgG_bam/"), pattern = "bai", invert =TRUE, value = TRUE))
sample.sheet$Control_ID <- rep(c("IgG"), times=4)
sample.sheet <- sample.sheet[,c(1,4,5,2,8,7,3,6)]

H3K36me3.dba <- dba(sampleSheet=sample.sheet)
H3K36me3.dba

H3K36me3.dba<- dba.count(H3K36me3.dba)
H3K36me3.dba


H3K36me3.dba <- dba.normalize(H3K36me3.dba)
norm <- dba.normalize(H3K36me3.dba, bRetrieve = TRUE)
norm

dba.plotPCA(H3K36me3.dba,  attributes=DBA_CONDITION, label=DBA_ID)

#Setting 37°C as the baseline condition
H3K36me3.dba <- dba.contrast(H3K36me3.dba, reorderMeta=list(Condition="37"), minMembers = 2)
H3K36me3.dba

H3K36me3.dba <- dba.analyze(H3K36me3.dba)
dba.show(H3K36me3.dba, bContrasts=TRUE)
#Pull up all peaks and plot on a VolcanoPlot
results <- dba.report(H3K36me3.dba, contrast = 1, th=1)
sum(results$FDR < .05)
sum(results$FDR < .1)
results
results.df <- as.data.frame(results)
write.csv(results.df, file = "NPC_H3K36me3_AllPeaks.csv", row.names=FALSE)
write.table(results.df, file = "NPC_H3K36me3_AllPeaks.bed", row.names=FALSE, sep="\t")

#Annotating Peaks
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

anno_function <- function(bedfile) {
  cuttag_anno <- annotatePeak(bedfile, TxDb = txdb, annoDb="org.Hs.eg.db")
  anno_df <- as.data.frame(cuttag_anno)
  plotAnnoPie(cuttag_anno)
  return(anno_df)
}

H3K36me3_results_anno <- anno_function(results)
write.csv(H3K36me3_results_anno, file = "NPC_H3K36me3_AllPeaks_anno.csv", row.names=FALSE)

#Plot all peaks on a volcanoplot
library("EnhancedVolcano")
library(dplyr)
#Plot p value vs fold change
EnhancedVolcano(H3K36me3_results_anno,
                lab = H3K36me3_results_anno$SYMBOL,
                x = 'Fold',
                y = 'p.value',
                xlim = c(-20, 20),
                labSize = 4.0,
                drawConnectors = TRUE, max.overlaps = Inf, widthConnectors = 0.25, arrowheads = TRUE,
                col=c("Black", "Dark grey", "Dark grey", "#A92003"), colAlpha = 1)
