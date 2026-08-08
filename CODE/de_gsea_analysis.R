### Analysis of ISCs with old mitochondria (Pop 2) vs ISCs with young mitochondria (Pop 1) ###

.libPaths("/scratch/project_2001351/Daniel/CONDA/ENVS/andersson1/lib/R/library")

library(tximport, quietly=T)
library(GenomicFeatures, quietly=T)
library(readr, quietly=T)
library(DESeq2, quietly=T)
library(RColorBrewer, quietly=T)
library(pheatmap, quietly=T)
library(limma, quietly=T)
library(biomaRt, quietly = T)
library(gplots, quietly = T)
library(ggplot2, quietly=T)
library(apeglm, quietly=T)

# Specify files

samples <- c("PK8207_Pop1","PK8207_Pop2","PK8208_Pop1","PK8208_Pop2","PK8209_Pop1","PK8209_Pop2","PK8293_Pop1","PK8293_Pop2","PK2323_ymito","PK2323_omito","PK2326_ymito","PK2326_omito")
files <- file.path("../DATA/COUNTS/",samples,"abundance.tsv")
names(files) <- samples
all(file.exists(files))

# Create transcript database

trdb <- makeTxDbFromGFF(file="../DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf", format="gtf")
k <- keys(trdb, keytype = "TXNAME")
tx2gene <- select(trdb, k, "GENEID", "TXNAME")

# tximport

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

saveRDS(txi.kallisto.tsv,file="../DATA/COUNTS/a3.txi.kallisto.rds")

##########
# DESeq2 #
##########

sampleTable <- data.frame(population = factor(rep(c("pop1", "pop2"),6)), mouse = factor(rep(c("PK8207","PK8208","PK8209","PK8293","PK2323","PK2326"), each = 2)))

rownames(sampleTable) <- colnames(txi.kallisto.tsv$counts)

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, ~mouse+population)

keep <- rowSums(counts(dds) >= 5) >= 6
table(keep)
dds <- dds[keep,]

dds <- DESeq(dds)

saveRDS(dds,"../DATA/COUNTS/a3.dds.rds")

############################
### Exploratory analysis ###
############################

vsd <- vst(dds, blind=FALSE)

# Sample distances

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$population, vsd$mouse, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
         
pdf("RESULTS/OLD_VS_YOUNG/sampleCorrelationHeatmap.pdf",6,5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# PCA

pcaData <- plotPCA(vsd, intgroup=c("population", "mouse"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=population, shape=mouse)) +
  geom_point(size=3) + scale_shape_manual(values = 1:nlevels(sampleTable$mouse)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
  
pdf("RESULTS/OLD_VS_YOUNG/pca.pdf",5,5)
ggplot(pcaData, aes(PC1, PC2, color=population, shape=mouse)) +
  geom_point(size=3) + scale_shape_manual(values = 1:nlevels(sampleTable$mouse)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd)[select,]))

pcaLoadings <- pca$rotation
head(sort(pcaLoadings[,1],decreasing=T),40)
head(sort(pcaLoadings[,2],decreasing=T),40)


# Cook's distances

boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

pdf("RESULTS/OLD_VS_YOUNG/cooks.pdf",5,5)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

# Dispersions

plotDispEsts(dds)

pdf("RESULTS/OLD_VS_YOUNG/dispersions.pdf",5,5)
plotDispEsts(dds)
dev.off()

### Gene complexity ###

library(SingleCellExperiment, quietly=T)
library(scater, quietly=T)

mitoGenes <- c("ENSMUSG00000064372.1","ENSMUSG00000064371.1","ENSMUSG00000064370.1","ENSMUSG00000064369.1","ENSMUSG00000064368.1","ENSMUSG00000064367.1","ENSMUSG00000064366.1","ENSMUSG00000064365.1","ENSMUSG00000064364.1",
"ENSMUSG00000064363.1","ENSMUSG00000065947.3","ENSMUSG00000064361.1","ENSMUSG00000064360.1","ENSMUSG00000064359.1","ENSMUSG00000064358.1","ENSMUSG00000064357.1","ENSMUSG00000064356.3","ENSMUSG00000064355.1",
"ENSMUSG00000064354.1","ENSMUSG00000064353.1","ENSMUSG00000064352.1","ENSMUSG00000064351.1","ENSMUSG00000064350.1","ENSMUSG00000064349.1","ENSMUSG00000064348.1","ENSMUSG00000064347.1","ENSMUSG00000064346.1",
"ENSMUSG00000064345.1","ENSMUSG00000064344.1","ENSMUSG00000064343.1","ENSMUSG00000064342.1","ENSMUSG00000064341.1","ENSMUSG00000064340.1","ENSMUSG00000064339.1","ENSMUSG00000064338.1","ENSMUSG00000064337.1",
"ENSMUSG00000064336.1")
mitoGenes <- mitoGenes[mitoGenes %in% rownames(dds)]
riboGenes <- read.table("../DATA/rRnaGenes.txt",header=T,stringsAsFactors=F,sep="\t")[,1]
riboGenes <- riboGenes[riboGenes %in% rownames(dds)]

sce <- SingleCellExperiment(assays=list(counts=counts(dds)))
sce$mouse <- dds$mouse
sce$sample <- colnames(sce)
sce$group <- dds$population

sce <- addPerCellQC(sce,subsets=list(mitochondrial=mitoGenes,rRNA=riboGenes))
sce$subsets_mitochondrial_percent
sce$subsets_rRNA_percent

# Convert ensembl IDs to gene names

geneIds <- rownames(dds)
G_list <- rtracklayer::import("../DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf")
G_list <- mcols(G_list[,c("gene_id","gene_name")])
colnames(G_list) <- c("ensembl_gene_id","external_gene_name")
G_list <- unique(G_list) 
gnames <- data.frame(ensembl_gene_id=geneIds)
gnames <- merge(gnames,G_list,by="ensembl_gene_id")
gnames$ensembl_gene_id <- as.character(gnames$ensembl_gene_id)
gnames$external_gene_name <- as.character(gnames$external_gene_name)
gnames$external_gene_name <- make.unique(gnames$external_gene_name) # Make duplicate gene names unique

# Add gene symbols to the sce object

rowData(sce)$symbol  <- rownames(dds)
rowData(sce)$symbol <- gnames$external_gene_name[match(rowData(sce)$symbol,gnames$ensembl_gene_id)]

# Plot cumulative gene contribution

pdf("RESULTS/OLD_VS_YOUNG/scaterPlot.pdf",10,5)
plotScater(sce,block1="group",colour_by = "mouse",nfeatures=1000)
dev.off()


########################################
### Differential expression analysis ###
########################################

res <- results(dds,alpha=0.1)

# Shrink log fold changes

resLFC <- lfcShrink(dds, coef="population_pop2_vs_pop1", type="apeglm")

# Add shrunken lfc estimates to results

res$lfcShrunken <- resLFC$log2FoldChange

# Convert ensembl gene ids to gene symbols and add them to the DE results table

res$symbol <- rownames(res)
res$symbol <- gnames$external_gene_name[match(res$symbol,gnames$ensembl_gene_id)]

# Visualize the p-value distribution

library(ggplot2, quietly=T)

ggplot(as.data.frame(res), aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0) + ggtitle("Wald test p-values")

# Sort by adjusted p-value

res <- res[order(res$padj),]

# Save the results table

write.table(res,"RESULTS/OLD_VS_YOUNG/degs.txt",quote = F,sep = "\t",row.names = T,col.names = NA)

# Visualize expression of DE genes

#res <- read.table("RESULTS/OLD_VS_YOUNG/degs.txt",header=T,stringsAsFactors=F,sep="\t",row.names=1)
degs <- rownames(res)[!is.na(res$padj) & res$padj <0.1]
degSymbols <- res$symbol[rownames(res) %in% degs]

for (i in 1:length(degs)) {

d <- plotCounts(dds, gene=degs[i], intgroup="population",main=degSymbols[i],returnData=TRUE)
d$mouse <- dds$mouse

pdf(paste0("RESULTS/OLD_VS_YOUNG/DE/",degSymbols[i],"_normcounts.pdf"),6,4)
print(ggplot(d, aes(x=population, y=count)) + geom_point(aes(shape=mouse, color=population),position=position_jitter(w=0.05,h=0),size=3) + ggtitle(degSymbols[i]) + geom_line(aes(group = mouse),size=0.1) + ylab("Normalized counts"))
dev.off()

}


############
### GSEA ###
############

#dds <- readRDS("../DATA/COUNTS/a3.dds.rds")

library(msigdbr, quietly = T)

# Load gene sets

gsets <- msigdbr(species = "Mus musculus") %>% dplyr::filter((gs_cat %in% c("H","C2")))

# Select genes seen in the data

gsets <- gsets[gsets$gene_symbol %in% gnames$external_gene_name,]

gsetlist <- gsets %>% split(x = .$gene_symbol, f = .$gs_name)

# Read IEC marker lists from Haber et al. 2017

eecMarkers <- read.table("../DATA/eecMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
enterocyteMarkers <- read.table("../DATA/enterocyteMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
epLateMarkers <- read.table("../DATA/epLateMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
gobletMarkers <- read.table("../DATA/gobletMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
panethMarkers <- read.table("../DATA/panethMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
stemMarkers <- read.table("../DATA/stemMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
tuftMarkers <- read.table("../DATA/tuftMarkers_haberEtAl2017_plate.txt",stringsAsFactors=F)$V1
mex3aMarkers <- read.table("../DATA/mex3aMarkers_barrigaEtAl2017.txt",stringsAsFactors=F)$V1
lrcSignature <- read.table("../DATA/lrcSignature_barrigaEtAl2017.txt",stringsAsFactors=F)$V1
ibpGenes <- read.table("../DATA/ibpGenes_kimEtAl2016.txt",stringsAsFactors=F)$V1

# Select markers seen in the data

mlist <- list(eecMarkers=eecMarkers,enterocyteMarkers=enterocyteMarkers,epLateMarkers=epLateMarkers,gobletMarkers=gobletMarkers,panethMarkers=panethMarkers,stemMarkers=stemMarkers,tuftMarkers=tuftMarkers,
mex3aMarkers=mex3aMarkers,lrcSignature=lrcSignature,ibpGenes=ibpGenes)
test <- lapply(mlist, FUN = function(x) {x <- x[x %in% gnames$external_gene_name]})

# Add marker lists to gsetlist

gsetlist <- c(gsetlist, mlist)

### Set gene indices ###

# Get voom-normalized values for genes

dE <- model.matrix(design(dds),colData(dds))
cE <- counts(dds,normalized=F)
voomE <- voom(cE,design=dE,plot=T) # check plot for the need to filter

# Check normalized expression value boxplots

boxplot(voomE$E)

# Convert ensembl gene ids to gene symbols

rnamesE <- rownames(voomE)
rnamesE <- gnames$external_gene_name[match(rnamesE,gnames$ensembl_gene_id)]

rownames(voomE) <- rnamesE

# Set indexes for GSEA

idxsE <- ids2indices(gsetlist,rownames(voomE),remove.empty = F)
names(idxsE) <- names(gsetlist)

# Remove sets smaller than 10 genes

idxsE <- idxsE[lengths(idxsE)>=10]

# Perform GSEA

set.seed(123)
gseaRes <- camera(voomE,idxsE,dE,use.ranks=T)

# Save results

write.table(gseaRes,"RESULTS/OLD_VS_YOUNG/GSEA/gsea_camera_results.txt",quote = F,sep = "\t",row.names = T,col.names = NA)


### Top candidates and visualization ###

#library(limma)
#gseaRes <- read.table("RESULTS/OLD_VS_YOUNG/GSEA/gsea_camera_results.txt",stringsAsFactors=F)
#res <- read.table("RESULTS/OLD_VS_YOUNG/degs.txt",stringsAsFactors=F)

# List top gene set candidates

upGsets <- c(rownames(gseaRes[gseaRes$FDR<0.1 & gseaRes$Direction=="Up",]),"panethMarkers")
downGsets <- rownames(gseaRes[gseaRes$FDR<0.1 & gseaRes$Direction=="Down",])

genesetsUp <- gsetlist[upGsets]
genesetsDown <- gsetlist[downGsets]

leadingUp <- genesetsUp

for (i in 1:length(genesetsUp)) {
  temp <- res[res$symbol %in% unlist(genesetsUp[i]),]
  leadingUp[i] <- list(head(temp$symbol[order(temp$stat,decreasing = T)],10))
  
}

leadingDown <- genesetsDown

for (i in 1:length(genesetsDown)) {
  temp <- res[res$symbol %in% unlist(genesetsDown[i]),]
  leadingDown[i] <- list(head(temp$symbol[order(temp$stat,decreasing = F)],10))
  
}

leading <- c(leadingUp,leadingDown)

{
sink("RESULTS/OLD_VS_YOUNG/GSEA/gsea_topCandidates.txt")
print(leading)
sink()
}

# Visualize

resIdxs <- ids2indices(gsetlist,res$symbol,remove.empty = F)
names(resIdxs) <- names(gsetlist)
resIdxs <- resIdxs[lengths(resIdxs)>0]

# Plot barcodeplots for select gene sets (stat on x-axis)

gsets <- c(genesetsDown,genesetsUp,mlist)

for (i in 1:length(gsets)) {
  pdf(paste0("RESULTS/OLD_VS_YOUNG/GSEA/geneset",i,"_stat.pdf"),5,4)
  barcodeplot(res$stat,index = unlist(resIdxs[names(gsets[i])]),main=strwrap(names(gsets[i])),cex.main=0.5)
  dev.off()
}
