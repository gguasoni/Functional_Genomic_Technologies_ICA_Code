#This code has been adapted from FGT_T7_Advanced_Analysis_RNASeq_UsingR


#The code for creating the boxplots was adapted from: https://r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html

#Load the required libraries & load the files for the tutorial
library(Rsubread)
library(edgeR)
library(limma)
library(gplots)
library(DESeq2)
library(affy)
library(QuasR)
library(pheatmap)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(scatterplot3d)


#Part 1: Quality Control
#This is the FastQC MultiQC stuff where you have to explain why you trimmed your data.

#Part2: Differential Gene Expression Analysis
#  This was done using the trimmed data. The "myfeaures_table" was made during the alignment by the professor. 

#This section is for loading and formatting all our base files

#Loading the annotation table that was prepared for us with the alignment
load("mytable_feaures")
#Create the gene x sample integer matrix
table_rnaseq<-(mytable_feaures)$counts
#Accessing column names from table_rnaseq
colnames_rnaseq <-colnames(table_rnaseq)
#Renaming the rows in the targets "adf" file to the SampleName
adf <-read.table("B188292_targets.csv",sep=',', row.names=1,fill=T,header=T)
#Ordering the elements in colnames_rnaseq to the rownames of adf
idx<-match(colnames_rnaseq,rownames(adf))
#Renaming the column names of table_rnaseq to the Sample name from adf for conciseness
colnames(table_rnaseq) <-adf$Sample[idx]
adf <-adf[idx,]
rownames(adf) <-colnames(table_rnaseq)
#viewing the re-formatted table_rnaseq
head(table_rnaseq)


#Statistical Analysis
#DESeq2 will used for the statistical analysis and the normalisation of the data. 

#This part is for loading the data into DESeq2 

#Creating the DESeq object for differential analysis
dds_rnaseq <- DESeqDataSetFromMatrix(countData = table_rnaseq, 
                                     colData = adf, 
                                     design = ~ treatment) #treatment is the only variable in this experiment
#Checking the object dimensions
dim(dds_rnaseq)
head(rownames(dds_rnaseq))
colnames(dds_rnaseq)
#check the object type
class(dds_rnaseq)
typeof(dds_rnaseq)

#Checking what this object contains
slotNames(dds_rnaseq)
colnames(colData(dds_rnaseq))
colnames(assay(dds_rnaseq))

#Quality Filtering and Normalisation of Data
Apply filters to remove genes that do not pass quality filters.


#first filter based upon counts
dds_rnaseq <- dds_rnaseq[rowSums(counts(dds_rnaseq) >= 10) >=3, ]
#We now have a DESeq object filtered for samples and genes that pass quality control 
nrow(dds_rnaseq)
#size factor estimation for normalisation- these are used to scale by library size when data is exported.
dds_rnaseq <- estimateSizeFactors(dds_rnaseq)
sizeFactors(dds_rnaseq)
#The size factors should relate to the library size.  We can check this from the annotation
adf_bases <- adf$Bases_.G.
names(adf_bases) <- rownames(adf)
adf_bases <- adf_bases/mean(adf_bases)
#Check how these two tables relate
plot(sizeFactors(dds_rnaseq), adf_bases)


## DESeq2 also provides methods that outputs data with various normalisations.
#Here we generate multiple different normalisations and compare between then- but for the rest of the tutorial we use rlog normalisation!  Note that the normalised data is not used for differential expression analysis- just for data visualisation.


#Normalised log counts
counts_rnaseq <- log2(1+counts(dds_rnaseq, normalized=TRUE))
#fpm is similar to cpm and discussed here https://rdrr.io/bioc/DESeq2/man/fpm.html.  Here the data is log transformed
fpm_rnaseq <-log(fpm(dds_rnaseq)+1)#this is also normalised by the sizeFactors
#rlog ad vsd are discussed here https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
vsd_rnaseq <- vst(dds_rnaseq, blind = T)
head(vsd_rnaseq, 3)
rld_rnaseq <- rlog(dds_rnaseq, blind = T)
head(rld_rnaseq, 3)



##Plots with Normalised Data.


adf_copy <- adf

par(mfrow = c(2, 2))

#custom_colours <- c("#D55E00", "#E69F00", "#F0E442", "#009E73","#0072B2", "#56B4E9", "#CC79A7")
custom_colours <- c("#D55E00", "#0072B2")

crd_as_df <- as.data.frame(counts_rnaseq)
crd_as_df <- rownames_to_column(crd_as_df, "Gene")
crd_as_df <- pivot_longer(crd_as_df, cols = -Gene, names_to = "Sample", values_to = "Count")
crd_as_df <- left_join(crd_as_df, adf_copy, by = "Sample")
crd_as_df$Colour <- custom_colours[crd_as_df$treatment]
crd_as_df$Colour <- custom_colours[as.numeric(as.factor(crd_as_df$treatment))]

boxplot(Count ~ Sample, data = crd_as_df,
        ylim = c(-2, 25),
        las = 2,          
        cex.axis = 0.7,
        xlab = "",
        ylab = "log2(Count)",
        main = "Boxplot: Count per Sample",
        col = crd_as_df$Colour)

fpm_as_df <- as.data.frame(fpm_rnaseq)
fpm_as_df <- rownames_to_column(fpm_as_df, "Gene")
fpm_as_df <- pivot_longer(fpm_as_df, cols = -Gene, names_to = "Sample", values_to = "Count")
fpm_as_df <- left_join(fpm_as_df, adf_copy, by = "Sample")
fpm_as_df$Colour <- custom_colours[fpm_as_df$treatment]
fpm_as_df$Colour <- custom_colours[as.numeric(as.factor(fpm_as_df$treatment))]

boxplot(Count ~ Sample, data = fpm_as_df,
        ylim = c(-2, 25),
        las = 2,           
        cex.axis = 0.7, 
        xlab = "",
        ylab = "FPM(Counts)",
        main = "Boxplot: FPM per Sample",
        col = fpm_as_df$Colour)

vsd_as_df <- as.data.frame(assay(vsd_rnaseq))
vsd_as_df <- rownames_to_column(vsd_as_df, "Gene")
vsd_as_df <- pivot_longer(vsd_as_df, cols = -Gene, names_to = "Sample", values_to = "Count")
vsd_as_df <- left_join(vsd_as_df, adf_copy, by = "Sample")
colour_levels <- unique(vsd_as_df$Sample)
vsd_as_df$Colour <- custom_colours[vsd_as_df$treatment]
vsd_as_df$Colour <- custom_colours[as.numeric(as.factor(vsd_as_df$treatment))]

boxplot(Count ~ Sample, data = vsd_as_df,
        ylim = c(-2, 25),
        las = 2,           
        cex.axis = 0.7,    
        main = "Boxplot: VSD per Sample",
        xlab = "", 
        ylab = "VST(Count)",
        col = vsd_as_df$Colour)

rld_as_df <- as.data.frame(assay(rld_rnaseq))
rld_as_df <- rownames_to_column(rld_as_df, "Gene")
rld_as_df <- pivot_longer(rld_as_df, cols = -Gene, names_to = "Sample", values_to = "Count")
rld_as_df <- left_join(rld_as_df, adf_copy, by = "Sample")
rld_as_df$Colour <- custom_colours[rld_as_df$treatment]
rld_as_df$Colour <- custom_colours[as.numeric(as.factor(rld_as_df$treatment))]

boxplot(Count ~ Sample, data = rld_as_df,
        ylim = c(-2, 25),
        las = 2,        
        cex.axis = 0.7,
        xlab = "",
        ylab = "RLD(Counts)",
        main = "Boxplot: RLD per Sample",
        col = rld_as_df$Colour)

#but suppose we look at non-normalised
na.rm=T
counts_un <-counts(dds_rnaseq, normalized=FALSE)
boxplot(log2(counts_un+1))
counts_un_rlog<-rlog(counts_un)
#We can make an mva.pairs plot for any of these data
#For example we can examine the unnormalized data.
#make mva.pairs-this is slow “mva.pairs.png” is provided 
#png(filename="mva.pairs.png",width=2000, height=2000)
mva.pairs(counts_un_rlog)
#dev.off()


#Exploratory Data Plots
#PCA plots similar to those we use for microarrays.


#we can run a PCA plot here...
#creating a Colour column
treatment <- colData(dds_rnaseq)$treatment
colour_vector <- ifelse(treatment == "Sham_operation", "#0072B2", "#D55E00")
#colour_vector <- ifelse(treatment == "Sham_operation", "#D55E00", "#0072B2")
colData(dds_rnaseq)$Colour <- colour_vector

legend_labels  <- c("Sham_operation", "BCAS")
legend_colours    <- c("#0072B2", "#D55E00")

# Perform PCA
pca <- prcomp(t(na.omit(assay(rld_rnaseq))), scale=T)
# Plot the PCA results
s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=colData(dds_rnaseq)$Colour)
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
#text(s3d.coords$x, s3d.coords$y, labels = colnames(rld_rnaseq),pos = 4,offset = 0.2,cex=0.5)
text(s3d.coords$x, s3d.coords$y, labels = colnames(rld_rnaseq),pos = sample(c(1,3), length(colnames(rld_rnaseq)), replace = TRUE),offset = 0.3,cex=0.5)

#or a quick plot in 2d
colour_palette <- c(Sham_operation = "#0072B2", `bilateral_carotid_artery stenosis_(BCAS)_operation` = "#D55E00")
qplot(pca$x[,1],pca$x[,2],
      xlab="PCA1", 
      ylab="PCA2",
      colour = treatment) + 
      geom_text(aes(label = colnames(rld_rnaseq)),      
                size = 2.2,
                vjust = -0.6, 
                show.legend = FALSE) +  
      scale_colour_manual(values = colour_palette,
                          name = "Treatment",
                          labels = c("BCAS", "Sham"))

#We can also output a distance matrix…
sampleDists <- dist(t(assay(rld_rnaseq)))
sampleDists

## Generate further heatmaps with the data

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix)
#easier PCA plot
plotPCA(rld_rnaseq, intgroup = c("treatment"))


## Differential Gene Expression
# Here we perform a two factor analysis using DESeq2.

dds_rnaseq <- DESeq(dds_rnaseq)
resultsNames(dds_rnaseq)
#This is how you specify a contrast…
result_treatment <-results(dds_rnaseq, name = "treatment_Sham_operation_vs_bilateral_carotid_artery.stenosis_.BCAS._operation")

result_treatment_df <- as.data.frame(result_treatment)
#dim(result_treatment_df) gives 18796
#to clean this df, we're going to get rid of rows where pval is NA, and rows where padj is NA. This is approx 370 genes, which is not a lot comparing to >18,000
result_treatment_df <- result_treatment_df[!is.na(result_treatment_df$pvalue),]
# dim(result_treatment_df) gives 18787
result_treatment_df <- result_treatment_df[!is.na(result_treatment_df$padj),]
#dim(result_treatment_df) gives 18422

#make a nice plot of differential expression
plotMA(result_treatment, main="DESeq2 dispersion plot", ylim=c(-2,2))
table(result_treatment$padj < 0.01)
table(result_treatment$padj < 0.05)
result_treatment_selected <- subset(result_treatment, padj < 0.05)
result_treatment_selected <-result_treatment_selected[order(abs(result_treatment_selected$log2FoldChange),decreasing = TRUE), ]
head(result_treatment_selected)

## Make heatmaps from each of the output tables- selected top 50 or all using fold change and FDR.

top10 <- rownames(result_treatment_selected)[1:10]
#make a heatmap for the top 50
pheatmap(assay(rld_rnaseq)[top10,],scale="row",show_rownames=T,main="RLD values, Row Scaled, Top 10")
#all
pheatmap(assay(rld_rnaseq)[rownames(result_treatment_selected),],scale="row",show_rownames=F,main="Treatment Selected, Row Scaled")



## Build Gene Annotation. A serialised result "mart" has been generated in case Ensembl is unavailable

#Download Ensembl annotation using BiomaRt and rename the samples
library(biomaRt)
#UK ensembl is being updated so we use a USA mirror, "useast.ensembl.org"
ensembl_host <-"https://www.ensembl.org"
head(biomaRt::listMarts(host = ensembl_host), 15)
head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL",host = ensembl_host))), 40) 
mart <- biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",mart    = useMart("ENSEMBL_MART_ENSEMBL",host    = ensembl_host))
resultAnnot <- biomaRt::getBM(values=rownames(dds_rnaseq),attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","description","strand"),filters="ensembl_gene_id",mart=mart)

#Merge Annotation with Input Data

#merge with input data
names <- resultAnnot[,1]
resultAnnot <- as.data.frame(resultAnnot)
rownames(resultAnnot) = names
idx<-match(rownames(dds_rnaseq),rownames(resultAnnot))
#make sure annotation is in same order
all(rownames(dds_rnaseq) == rownames(resultAnnot))
grr<-resultAnnot[match(rownames(dds_rnaseq), resultAnnot$ensembl_gene_id),]
all(rownames(dds_rnaseq) == rownames(grr))
resultAnnot <-grr
all(rownames(dds_rnaseq) == rownames(resultAnnot))
#make the nice names
nice_names<-paste(resultAnnot$ensembl_gene_id,resultAnnot$external_gene_name, sep = '_')
resultAnnot$nice_names <-nice_names
head(resultAnnot)
all(rownames(dds_rnaseq) == rownames(resultAnnot))
#check names
rld_rnaseq <- rlog(dds_rnaseq, blind = TRUE)#
idx2 <-match(rownames(result_treatment_selected)[1:10],rownames(dds_rnaseq))
plotme <-(rld_rnaseq)[rownames(result_treatment_selected)[1:10],]
rownames(plotme)<-resultAnnot$nice_names[idx2]
#make heatmap with candidate genes
png(filename="heatmap_candidates.png",width=800, height=1000)
pheatmap(assay(plotme),scale="row",fontsize_row = 10,cellheight =12, cellwidth =12,treeheight_row = 40, treeheight_col = 40)
dev.off()

pheatmap(assay(plotme),scale="row",fontsize_row = 10,cellheight =12, cellwidth =12,treeheight_row = 40, treeheight_col = 40)

#dev.off()


## Functional Enrichment Analysis


#Downloading the fgsea table because we need it now
library(fgsea)

download.file("https://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata", destfile = "mouse_H_v5p2.rdata")
load("mouse_H_v5p2.rdata")

#We need to map to EntrezIDs- the annotation uses EntrezID
head(Mm.H)

#
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
#listAttributes(mart)
ens2entrez <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), mart=mart)
head(ens2entrez)

#Map gene names to entrezgene_id
#if NA remove (we only need to consider mappable genes)
result_treatment_df$names2 <- rownames(result_treatment_df)

results_treatment2 <- inner_join(result_treatment_df, ens2entrez, by=join_by("names2"=="ensembl_gene_id"))

#filtering for genes where entrez_gene_id = NA
results_treatment2 <- results_treatment2[!is.na(results_treatment2$entrezgene_id),]

results_treatment2 <-results_treatment2[order(results_treatment2$pvalue,decreasing=TRUE),]

results_treatment3 <-results_treatment2[,c("pvalue","entrezgene_id")]

#remove non-unique row names- we don't mind deleting these in this case
#it is wise to check the genes that are duplicated to see if they should be excluded
results_treatment3 <-results_treatment3[!duplicated(results_treatment3$entrezgene_id),]

#results_treatment4 <- as_tibble(results_treatment3)
#results_treatment4 <- column_to_rownames(results_treatment4, var = "entrezgene_id")

rownames(results_treatment3)<-results_treatment3[,"entrezgene_id"]
results_treatment3["entrezgene_id"] <-NULL

#it might want  a named vector
original_entrez_ids <-rownames(results_treatment3)
colnames(results_treatment3) <-NULL
results_treatment4 <-results_treatment3[,1]
names(results_treatment4)<-original_entrez_ids

#time to perform the enrichment analysis
results_treatment5  <- results_treatment4[!(is.na(names(results_treatment4)))]
results_treatment5  <- results_treatment5[!(is.na(results_treatment5))]
results_treatment5 <-log2(results_treatment5)

fgseaRes <- fgsea(Mm.H, results_treatment5, minSize=25, maxSize = 500)

fgseaRes <-fgseaRes[order(fgseaRes$NES,decreasing=FALSE),]
head(fgseaRes)

plotEnrichment(Mm.H[["HALLMARK_GLYCOLYSIS"]],results_treatment5)

#plotting the enrichment data
top10_fgsea_res <- fgseaRes[fgseaRes$padj < 0.05, ]
top10_fgsea_res <- top10_fgsea_res[order(top10_fgsea_res$NES, decreasing = TRUE), ]
top10_fgsea_res <- head(top10_fgsea_res, 10)
top10_fgsea_res$pathway <- factor(top10_fgsea_res$pathway, levels = top10_fgsea_res$pathway[order(top10_fgsea_res$NES)])

ggplot(top10_fgsea_res, aes(x=NES, y=pathway)) +
  geom_col(aes(fill=top10_fgsea_res$padj)) +
  scale_fill_viridis_c(option = "B") +
  theme_minimal() +
  labs(title="Top Enriched Pathways", x="Normalized Enrichment Score", y="Pathway")


