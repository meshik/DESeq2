#####DESeq2 analysis##### (requires at least 2 biological replicates per condition)
# legend: # <- a note on what is going on
          ## <- describe the next few lines of code
          ### <- abandoned code, left to rot forever
# make sure you get no errors when loading libraries. If you do, then find an install command via Google
library("tximport") #
library("readr")
library(DESeq2)
library(biomaRt)
library(stringr)
library("IHW")
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(gridExtra)
library(grid)
library("RColorBrewer")
library(dplyr)
library(fgsea)
library(goseq)  # for GO enrichment analysis
library(matrixStats)  # for rowMedians()
library(bcbioRNASeq)  # for plotCorrelationHeatmap()


## Specify the full paths to all of your ".genes.results" files (from RSEM), and save those paths in 'files'
dir <- '/Users/You/Desktop/Analysis/rsem results/'  # the folder with your RSEM output
samples <- c("sample1", "sample2", "sample3")  # the names of your samples
files <- file.path(dir2, samples2, ".genes.results", fsep = '')
names(files) <- c("24h_condition", "24h_condition", "4h_condition", "4h_condition", "Control")  # name all your conditions in the same order as two lines above
all(file.exists(files))  # Use this line to see if your path is OK: prints TRUE if all files do exist; FALSE otherwise

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)  # produces an initial count table
head(txi.rsem$counts)  # head() prints the first few rows of the table


## Prepare the DESeq2 object
sampleTable <- data.frame(condition = factor(c("24h_condition", "24h_condition", "4h_condition", "4h_condition", "Control")))
txi.rsem$length[txi.rsem$length == 0] <- 1  # convert all counts with 0 to be 1
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

keep <- rowSums(counts(dds)) >= 10  # pre-filtering: anything with less than 10 counts is thrown away
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("24h_condition", "4h_condition", "Control"))  # enter each condition category once here
dds <- DESeq(dds)  # produces the DESeqDataSet object


## biomaRt - convert gene IDs from ENSEMBL to symbol or Entrez. Code is written for mouse
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")  # "ENSEMBL_MART_ENSEMBL" chosen from listMarts()
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)  # "mmusculus_gene_ensembl" chosen from listDatasets(ensembl)

filterType <- "ensembl_gene_id"
filterValues <- rownames(dds2)  # the gene IDs we want to change
filterValues2 = sapply(strsplit(filterValues, '.', fixed=T), function(x) x[1])  # this removes the version extension of the ensembl gene ID

attributes = listAttributes(ensembl) 
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol')

annot <- getBM(attributes = attributeNames, filters = filterType, 
               values = filterValues2, mart = ensembl, uniqueRows = TRUE)

isDup <- duplicated(annot$ensembl_gene_id)  # boolean (True/False) vector of duplicates
dup <- annot$ensembl_gene_id[isDup]  # ENSEMBL gene IDs of the duplicated IDs
annot <- annot[!isDup,]  # this removes all of the duplicate mappings caused by converting gene IDs

annot$entrezgene_id[is.na(annot$entrezgene_id)] <- "999999999"  # changes all of the NA to 999999999
annot$mgi_symbol[is.na(annot$mgi_symbol)] <- "notRealGene"  # change all of the NA to "notRealGene"
# if you see in your data "notRealGene" or 999999999, then these are genes you SHOULD NOT include in any analysis

rownames(dds) <- annot$mgi_symbol  # change the gene IDs in our DESeqDataSet to symbol. You can type in "eentrezgene_id" instead of "mgi_symbol"


## Create the condition contrasts. You should create a "res" (DESeqResults) object for each comparison you want to observe
res24v0 <- results(dds, contrast=c("condition", "24h_condition", "Control"))
res24v4 <- results(dds, contrast=c("condition", "24h_condition", "4h_condition"))
res4v0 <- results(dds, contrast=c("condition", "4h_condition", "Control"))

# corrects res for log-fold-change
### resLFC2 <- lfcShrink(dds2, coef="condition_Control_vs_24h_Irisin", type="apeglm")

# how many adjusted p-values are under 0.05(or anything else)?:
sum(res24v0$padj < 1e-10, res24v4$padj < 1e-10, res4v0$padj < 1e-10, na.rm=TRUE)

# also does some normalization
resIHW2 <- results(dds2, filterFun=ihw)
summary(resIHW2)

## PLOT: MA
par(mfrow=c(3,1)) # creates a m*n plot array for the following plots
plotMA(res24v0, main = "24h vs Ctrl"); plotMA(res24v4, main = "24h vs 4h"); plotMA(res4v0, main = "4h vs Ctrl")

# can use this next command to click on the plot and select points of interest. Press 'escape' to exit this mode
# idx2 <- identify(res2$baseMean, res2$log2FoldChange)
# rownames(res2)[idx2]

## PLOT: compare single genes across the conditions
plotCounts(dds2, gene="Plxdc2", intgroup="condition")

# this gives log2(n + 1)
ntd2 <- normTransform(dds2)


## FUNCTION: takes a list of genes and automatically copies it to your clipboard, each on it's own line
copy_to_clipboard = function(x,sep="\t",col.names=T,...) { 
  write.table(x
              ,file = pipe("pbcopy")
              ,sep=sep
              ,col.names = col.names
              ,row.names = F
              ,quote = F,...)
}
# copy_to_clipboard(deGenes1)


## PLOT: pheatmap
# define a threshold and take all the gene names (@rownames) that surpass that threshold
#TODO: make a function that takes any amount of DESeqResults and merges their differentially expressed genes into one vector...
# first need to understand how to enter any amount of objects as input arguments...
# also need to for loop over each pair and do the union of them
deThreshold <- 1e-15
deGenes1 <- res24v0@rownames[res24v0$padj < deThreshold]
deGenes1 <- deGenes1[!is.na(deGenes1)]  # good command to know to remove all na from a vector
deGenes2 <- res24v4@rownames[res24v4$padj < deThreshold]
deGenes2 <- deGenes2[!is.na(deGenes2)]
deGenes3 <- res4v0@rownames[res4v0$padj < deThreshold]
deGenes3 <- deGenes3[!is.na(deGenes3)]

# take the union of all the DE genes
uni <- union(deGenes1, deGenes2) 
uni <- union(uni, deGenes3)
uni <- uni[uni != "999999999" & uni != ""]

## PLOT: pheatmap - this function has many possible arguments that are all useful and easy to use
cdata2 <- colData(dds2)
phm <- pheatmap(assay(ntd2)[uni,], cluster_rows=TRUE, show_rownames=FALSE,
                cluster_cols=FALSE, scale = "row", clustering_distance_rows = "euclidean",
                main = "Union of differentially expressed genes (log2(n+1)) from each pair of conditions.\n p-adj threshold: 0.001\n# of genes: 1201",
                annotation_col=as.data.frame(cdata2[,"condition"], row.names=rownames(cdata2)), cutree_rows = 5)
### arguments for pheatmap: cutree_rows = 7, scale = "row",

# This next section is dependent on first running Seurat::FindMarkers() on a Seurat object
# mghm <- pheatmap(assay(ntd2)[rownames(mgDiffExp[c(1:5, 7:36, 38:40),]),], cluster_rows=TRUE, show_rownames=TRUE,
#                 cluster_cols=FALSE, scale = "row", main = "Migroglia DE genes against all",
#                 annotation_col=as.data.frame(cdata2[,"condition"], row.names=rownames(cdata2)))
# astrohm <- pheatmap(assay(ntd2)[rownames(astroDiffExp[c(1:8, 10:21, 23:31, 33:40),]),], cluster_rows=TRUE, show_rownames=TRUE,
#                     cluster_cols=FALSE, scale = "row", main = "Astrocyte DE genes against all",
#                     annotation_col=as.data.frame(cdata2[,"condition"], row.names=rownames(cdata2)))
# avmhm <- pheatmap(assay(ntd2)[rownames(bothDiffExpAvM[c(1:19, 21:25, 27:40),]),], cluster_rows=TRUE, show_rownames=TRUE,
#                    cluster_cols=FALSE, scale = "row", main = "Astrocyte vs Microglia DE genes",
#                    annotation_col=as.data.frame(cdata2[,"condition"], row.names=rownames(cdata2)))
# mvahm <- pheatmap(assay(ntd2)[rownames(bothDiffExpMvA[c(1:9, 11:26, 28:40),]),], cluster_rows=TRUE, show_rownames=TRUE,
#                   cluster_cols=FALSE, scale = "row", main = "Microglia vs Astrocyte DE genes",
#                   annotation_col=as.data.frame(cdata2[,"condition"], row.names=rownames(cdata2)))

# Extract data from pheatmap row clustering in order to inspect gene groups (modules)
plot(phm$tree_row)
abline(h=2, col="red", lty=2, lwd=2)

clustGenes <- sort(cutree(phm$tree_row, h=4.5))  # cuts the rows into clusters, dependent on the height specified (h)

# copy to clipboard any cluster you want to examine
copy_to_clipboard(names(clustGenes[clustGenes==1]))


## PLOT: PCA
vsd <- vst(dds2, blind=FALSE)  # variance stabilizing transformation by closed-form expression
plotPCA(vsd, intgroup="condition")

## PLOT: Pearson correlation heatmap between samples
plotCorrelationHeatmap(vsd, method = "pearson")

## PLOT: Sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,col=colors)

## PLOT: volcano
# create individual plot and save them
volc1 <- EnhancedVolcano(res24v0, lab = rownames(res24v0), x = 'log2FoldChange',
                         y = 'pvalue', xlim = c(-3, 3), ylim = c(0, 75), 
                         title = "24h vs Control", legendPosition = 'bottom') 
volc2 <- EnhancedVolcano(res24v4, lab = rownames(res24v4), x = 'log2FoldChange',
                         y = 'pvalue', xlim = c(-3, 3), ylim = c(0, 75), 
                         title = "24h vs 4h", legendPosition = 'bottom')
volc3 <- EnhancedVolcano(res4v0, lab = rownames(res4v0), x = 'log2FoldChange',
                         y = 'pvalue', xlim = c(-3, 3), ylim = c(0, 75), 
                         title = "4h vs Control", legendPosition = 'bottom')
# now we'll add all these plots to one figure
grid.arrange(volc1, volc2, volc3, ncol=3,
             top = textGrob('EnhancedVolcano', just = c('center'), gp = gpar(fontsize = 32)))

# Extract the significant genes on the volcano plot
signifVolcano <- function(res, pval=1e-5, lfc=1){
  return(rownames(subset(res, -log(res$padj) > 5 & (res$log2FoldChange > 1 | res$log2FoldChange < -1))))
}
volcGenes24v0 <- signifVolcano(res24v0)


## Enrichment Analysis using clusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db)
clusterProfiler <- function(geneSet, padjThresh=1e-3){
  # First take only DE genes from tghe DESeqResults input:
  geneSet <- na.omit(geneSet)  # remove all rows that have NA in any column
  geneSet <- geneSet[geneSet$padj < padjThresh,]  # take only DE rows
  geneSet <- geneSet[geneSet@rownames != '999999999',]  # removes the conversion to ENTREZ(NCBI) IDs
  
  geneList <- as.vector(geneSet$log2FoldChange)
  names(geneList) <- geneSet@rownames
  
  ggo <- clusterProfiler::groupGO(gene = names(geneList), OrgDb = org.Mm.eg.db, ont = "BP", level = 3, readable = TRUE)
  print(barplot(ggo, drop=TRUE, showCategory=24))  # inside functions, you need to explicitly print() any plot that isn't the last plot
  
  ego <- clusterProfiler::enrichGO(gene = names(geneList), OrgDb = org.Mm.eg.db, ont = "BP",
                                   pAdjustMethod = "BH",pvalueCutoff = 1e-4, qvalueCutoff = 1e-4, readable = TRUE)
  print(clusterProfiler::dotplot(ego, showCategory=25))
  
  kk1 <- clusterProfiler::enrichKEGG(gene = names(geneList), organism = 'mmu', pAdjustMethod = "BH",
                                     pvalueCutoff = 1e-5, qvalueCutoff  = 1e-5)
  print(cnetplot(kk1, categorySize="geneNum", foldChange=geneList))
  
  paste0('Number of genes < ', padjThresh, ' is: ', length(geneSet@rownames))
}

clusterProfiler(geneSet = res24v0); clusterProfiler(geneSet = res24v4); clusterProfiler(geneSet = res4v0)


## Gene Set Enrichment Analysis
GSEA <- function(res){
  ranks = res$log2FoldChange
  names(ranks) = res@rownames
  ranks = sort(ranks, decreasing = TRUE)
  ranks = ranks[!duplicated(names(ranks))]
  ranks = ranks[names(ranks) != 999999999]
  print(head(ranks))
  
  load("/Users/dekel/Desktop/020820/mouse_H_v5p2.rdata")  # saved as the variable 'Mm.H'
  fgseaRes <- fgsea(Mm.H, ranks, minSize=15, maxSize = 500, eps = 0)
  print(head(fgseaRes[order(padj, -abs(NES)), ], n=10))
  
  print(plot.new())  # added this because the plots tend to overlap
  
  topUp <- fgseaRes %>% 
    filter(ES > 0) %>% 
    top_n(10, wt=-padj)
  topDown <- fgseaRes %>% 
    filter(ES < 0) %>% 
    top_n(10, wt=-padj)
  topPathways <- bind_rows(topUp, topDown) %>% 
    arrange(-ES)
  
  print(plotGseaTable(Mm.H[topPathways$pathway], ranks, fgseaRes, gseaParam = 0.5))
  
  isSigGene <- res$padj < 0.01 & !is.na(res$padj) & res@rownames != '999999999'
  genes <- as.integer(isSigGene)
  nonDup <- res@rownames[res@rownames != '999999999']
  names(genes) <- nonDup
  print(head(nonDup))
  print(head(genes))
  
  pwf <- nullp(genes, "mm10", "ensGene", bias.data = rowMedians(dds2@assays@data$avgTxLength))
}

GSEA(res24v0); GSEA(res24v4); GSEA(res4v0)  # makes a plot of pathway enrichment of up and down regulated pathways


## Export DESeqResults to CSV files
col24 <- rowMeans(vsd@assays@data@listData[[1]][,1:3])
col4 <- rowMeans(vsd@assays@data@listData[[1]][,4:5])
col0 <- rowMeans(vsd@assays@data@listData[[1]][,6:8])
allCol <- cbind(col0, col4, col24)
write.csv(as.data.frame(allCol), file="timeCourse.csv")
# write.csv(as.data.frame(res24v0), file="res24v0.csv")
# write.csv(as.data.frame(res24v4), file="res24v4.csv")
# write.csv(as.data.frame(res4v0), file="res4v0.csv")


## Markers from literature:
damMarkers <- c("Cst7","Lpl","Clec7a","Itgax","Spp1","Igf1","Apoe","Axl",
                "Ank","Ch25h","Ctsd","Ccl3","Baiap2l2","Csf1","Fam20c","Gnas",
                "Ctsb","Lyz2","Cd63","Lgi2","Tyrobp","Hpse","P2ry12",
                "Mamdc2","Cd63-ps","Ccl6","Gm1673","Cd9","Ccl4","Ctsz",
                "Tmem119","Fth1","B2m","Cx3cr1","Gusb","Fabp3","Actr3b","H2-D1",
                "Ctsl","Serpine2","Selplg","Cd52","Lgals3bp","Dkk2",
                "Psat1","Gm10076","Cxcl14","Lox","Flt1","Capg","Mif","Nceh1",
                "Sulf2","Serinc3","Lyz1","Rps2","Npc2","Gm15500",
                "Etl4","Cd68","Rps28","Syngr1","Rpsa-ps10","Anxa5","Hif1a",
                "H2-K1","Atp1a3","Kcnj2","Hexa","Gpnmb","Cd34","Fabp5")
daaDamMarkers <- c("Prdx1","Ftl1","Npc2","Apoe","Aplp2","Scd2","Mt1","Fth1","Ctsd","Ctsb","H2-D1","B2m","Cd63","Ctsl","Lgals3bp","H2-K1","Cd9","Sulf2")
a1Markers <- c("H2-D1","Serping1","H2-T23","Ggta1","Iigp1","Gbp2","Fkbp5","Psmb8","Srgn","Amigo2")
a2Markers <- c("Clcf1","Ptx3","S100a10","Sphk1","Cd109","Ptgs2","Emp1","Tm4sf1","B3gnt5","Cd14")
panMarkers <- c("Lcn2","Steap4","S1pr3","Timp1","Hspb1","Cxcl10","Osmr","Cp","Aspg","Gfap","Vim")
phm <- pheatmap(assay(ntd2)[c("Cdk4", "Cdk6", "Cdk2", "Ccna2", "Ccnb1", "Ccnd1", "Cdc25a", "Cdc25b", "Cdc25c", "Cdt1"),],
                cluster_rows=FALSE, show_rownames=TRUE,
                cluster_cols=FALSE, clustering_distance_rows = "euclidean", scale = "row",
                annotation_col=as.data.frame(cdata2[,"condition"], row.names=rownames(cdata2)))
# , scale = "row"
ntd2[panMarkers]@assays@data@listData