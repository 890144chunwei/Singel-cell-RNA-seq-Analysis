install.packages("Seurat")
install.packages('Matirx')
install.packages("Matrix")
install.packages("tidyverse")
BiocManager::install('SingleCellExperiment')
BiocManager::install(c('scuttle', 'scran', 'scater', 'uwot', 'rtracklayer'))
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("limma", force = TRUE)
BiocManager::install("edgeR", force = TRUE)
BiocManager::install("scRNAseq")
library('Seurat')
library('dplyr')
library('ggplot2')
library("SingleCellExperiment")
library("scuttle")
library("scran")
library("scater")
library("tidyverse")
library("SingleR")
library("celldex")
library("org.Mm.eg.db")
library("edgeR")
library("scRNAseq")
library("pheatmap")
library(enrichplot)
library(clusterProfiler)

#Create matrix using SingleCellExperiment function
matrix_dir = "~/Desktop/scRNA-seq/outs/count/filtered_feature_bc_matrix/"
sce <- read10xCounts(matrix_dir)
assay(sce,"counts")
#gene annotation with mm10
path<-tempfile()
bfc <- BiocFileCache(path, ask = FALSE)
mm10.gtf <- bfcrpath(bfc, file.path("http://ftp.ensembl.org/pub/release-98","gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz"))
gene.data<- rtracklayer::import(mm10.gtf)
gene.data<- gene.data[gene.data$type=="gene"]
names(gene.data)<- gene.data$gene_id
is.gene.related <- grep("gene_",colnames(mcols(gene.data)))
mcols(gene.data) <- mcols(gene.data)[,is.gene.related]
rowRanges(sce)<- gene.data[rownames(sce)]
rowRanges(sce)[1:10,]
#Label samples from CellRanger aggr
aaa <- gsub(pattern= ".*-1", replacement="WT1", x= sce$Barcode)
aaa <- gsub(pattern= ".*-2", replacement="WT2", x= aaa)
aaa <- gsub(pattern= ".*-3", replacement="Mut1", x= aaa)
aaa <- gsub(pattern= ".*-4", replacement="Mut2", x= aaa)
sce$Genotype <-aaa
rm(aaa)

#basic  QC filtering mito-high
sce.unfilter <- sce
location<-rowRanges(sce.unfilter)
is.mito <- which(seqnames(location)=="MT")
df <- perCellQCMetrics(sce.unfilter,subsets=list(Mito=is.mito))
summary(df$sum)
high.mito <- isOutlier(df$subsets_Mito_percent, type = "higher", log = TRUE)
summary(high.mito)
summary(df$subsets_Mito_percent)
sce.unfilter$Mito.percent <- df$subsets_Mito_percent
##QC with fixed threshold
qc.feature <- df$sum > 30000
qc.count <- df$detected > 6500
qc.mito <- df$subsets_Mito_percent >10
discard <- qc.feature | qc.count | qc.mito
DataFrame(nFeature=sum(qc.feature), nCount=sum(qc.count), percentMito=sum(qc.mito), Total=sum(discard))
sce$discard <- discard
#Plotting QC metrics
gridExtra::grid.arrange(
plotColData(sce, x="Genotype", y="sum",colour_by = "discard")+ scale_y_log10()+ ggtitle("Total count"),
plotColData(sce, x="Genotype", y="detected",colour_by = "discard")+ ggtitle("Detected features"),
plotColData(sce, x="Genotype", y="subsets_Mito_percent",colour_by = "discard")+ ggtitle("Percent mito"),
ncol=1)
#Subsetting based on QC - Filter low-quality cells
sce.filter <- sce.unfilter[,df$subsets_Mito_percent <10]
sce.filter <- sce[,!discard]
#Normalization of library size factor
lib.sce.filter <- librarySizeFactors(sce.filter)
hist(log10(lib.sce.filter), xlab="Log10[Size Factor]", col='grey80')
#Scaling and log-transformation
sce.filter <- logNormCounts(sce.filter)
#Feature selection - Quantify per-gene variation
dec.sce.filter <- modelGeneVar(sce.filter)
fit.sce.filter <- metadata(dec.sce.filter)
plot(fit.sce.filter$mean, fit.sce.filter$var, pch=16, xlab="Mean of log expression", ylab="Variance of log expression")
curve(fit.sce.filter$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dec.sce.filter[order(dec.sce.filter$bio, decreasing = TRUE),]
#Selecting highly variable genes for PCA
hvg.sce.filter <- getTopHVGs(dec.sce.filter, prop = 0.2)

#Dimensional reduction - PCA
set.seed(100)
sce.filter <- fixedPCA(sce.filter, subset.row = hvg.sce.filter)
percent.var <- attr(reducedDim(sce.filter), "percentVar")
chosen.elbow <- findElbowPoint(percent.var)
plot(percent.var, xlab='PC', ylab='Variance explained (%)')
abline(v= chosen.elbow, col="red")
reducedDim(sce.filter, "PCA.elbow") <- reducedDim(sce.filter)[,1:chosen.elbow]
gridExtra::grid.arrange(
  plotReducedDim(sce.filter, dimred="PCA.elbow", colour_by = "Genotype"),
  plotReducedDim(sce.filter, dimred="PCA.elbow", colour_by = "Mito.percent"),
  ncol=2)
  
#Dimensional reduction - tSNE
set.seed(00101001101)
sce.filter<- runTSNE(sce.filter, dimred="PCA.elbow", perplexity=60)
#Dimensional reduction - UMAP
sce.filter<- runUMAP(sce.filter, dimred="PCA.elbow", n_neighbors=4)

#clustering and adjusting parameters
set.seed(100)
cluster.kmeans <- clusterCells(sce.filter, use.dimred = "PCA.elbow", BLUSPARAM = MbkmeansParam(center=13))
table(cluster.kmeans)
colLabels(sce.filter) <- cluster.kmeans

#Check dimensional reduction and clustering parameters
plotReducedDim(sce.filter, dimred="UMAP", colour_by = "label",  text_by = "label", text_colour = "red")
row.names(sce.filter) <- rowData(sce.filter)$gene_name
plotExpression(sce.filter, features= c("Ly6a","Kit"), x="label", colour_by = "label")

#Select reference profile with SingleR
mouse.immGen <- ImmGenData()
mouse.rnaSeq <- MouseRNAseqData()
counts <- as.matrix(assay(sce.filter,"counts"))
common <- intersect(rownames(counts), rownames(mouse.immGen))
common <- intersect(common, rownames(mouse.rnaSeq))
mouse.rnaSeq <- mouse.rnaSeq[common,]
mouse.immGen <- mouse.immGen[common,]
sce.filter.common <- sce.filter[common,]
sce.filter.common <- logNormCounts(sce.filter.common)
pred.immgen.main <- SingleR(test = sce.filter.common, ref = mouse.immGen, labels = mouse.immGen$label.main, assay.type.test = 1)
sort(table(pred.immgen.main$labels), decreasing= TRUE)
sce.filter$pred.immgen.main <- pred.immgen.main

#Subsetting
dump <- which(pred.immgen.main$labels =="Stem cells")
sce.filter.stem <- sce.filter[,dump]
lib.sce.filter.stem<- librarySizeFactors(sce.filter.stem)
sce.filter.stem <- logNormCounts(sce.filter.stem)
dec.sce.filter.stem <- modelGeneVar(sce.filter.stem)
fit.sce.filter.stem <- metadata(dec.sce.filter.stem)
hvg.sce.filter.stem <- getTopHVGs(dec.sce.filter.stem, prop = 0.2)
set.seed(100)
sce.filter.stem <- fixedPCA(sce.filter.stem, subset.row = hvg.sce.filter.stem)
percent.var <- attr(reducedDim(sce.filter.stem), "percentVar")
chosen.elbow <- findElbowPoint(percent.var)
reducedDim(sce.filter.stem, "PCA.elbow") <- reducedDim(sce.filter.stem)[,1:chosen.elbow]
sce.filter.stem<- runUMAP(sce.filter.stem, dimred="PCA.elbow", n_neighbors=4)
set.seed(100)
cluster.kmeans <- clusterCells(sce.filter.stem, use.dimred = "PCA.elbow", full=TRUE)
table(cluster.kmeans$clusters)
colLabels(sce.filter.stem) <- cluster.kmeans$clusters
dump <- which(sce.filter.stem$label != "20" & sce.filter.stem$label != "21" & sce.filter.stem$label != "22")
sce.filter.stem <- sce.filter.stem[,dump]
cluster.kmeans <- clusterCells(sce.filter.stem, use.dimred = "PCA.elbow", full=TRUE)
table(cluster.kmeans$clusters)
colLabels(sce.filter.stem) <- cluster.kmeans$clusters

#Visualize clustering after subseting
gridExtra::grid.arrange(
    plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "label",  text_by = "label", text_colour = "black"),
    plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Flt3",  text_by = "label", text_colour = "black"),
    plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Ly6a",  text_by = "label", text_colour = "black"),
    plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Gata2",  text_by = "label", text_colour = "black"),
  ncol=2)
sce.filter.common <- sce.filter.stem[common,]
pred.immgen.fine <- SingleR(test = sce.filter.common, ref = mouse.immGen, labels = mouse.immGen$label.fine, assay.type.test = 1)
sort(table(pred.immgen.fine$labels), decreasing= TRUE)
sce.filter.stem$pred.immgen.fine <- pred.immgen.fine
#clustering optimization
sil.approx <- approxSilhouette(reducedDim(sce.filter.stem,"PCA.elbow"), clusters = colLabels(sce.filter.stem))
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width>0, colLabels(sce.filter.stem), sil.data$other))
sil.data$cluster <- colLabels(sce.filter.stem)
ggplot(sil.data, aes(x=cluster, y=width, colour= closest))+ ggbeeswarm::geom_quasirandom(method="smiley")
table(Cluster=colLabels(sce.filter.stem), sil.data$closest)
pure.stem <- neighborPurity(reducedDim(sce.filter.stem,"PCA.elbow"),colLabels(sce.filter.stem))
pure.data <- as.data.frame(pure.stem)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- colLabels(sce.filter.stem)
ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +ggbeeswarm::geom_quasirandom(method="smiley")
table(Cluster=colLabels(sce.filter.stem), pure.data$maximum)

#Modularity
g <- cluster.kmeans$objects$graph
ratio<- pairwiseModularity(g, colLabels(sce.filter.stem), as.ratio = TRUE)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE, color=colorRampPalette(c("white", "blue"))(100))
cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1), mode="upper", weighted=TRUE, diag=FALSE)
set.seed(11001010)
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*5,layout=igraph::layout_with_lgl)
tab <- table(pred.immgen.fine$pruned.labels, sce.filter.stem$label)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
#ARI
cluster.5 <- clusterCells(sce.filter.stem, use.dimred = "PCA.elbow", BLUSPARAM=NNGraphParam(k=5))
cluster.10 <- clusterCells(sce.filter.stem, use.dimred = "PCA.elbow", BLUSPARAM=NNGraphParam(k=10))
pairwiseRand(cluster.10, cluster.5, mode="index")
breakdown <- pairwiseRand(ref=cluster.10, alt=cluster.5, mode="ratio")
pheatmap(breakdown, color=viridis::magma(100), cluster_rows=FALSE, cluster_cols=FALSE)
#Clustering annotation
plotExpression(sce.filter.stem, x="label", features="Ly6a")
gridExtra::grid.arrange(
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "label",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Slamf1",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Ly6a",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Flt3",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Spi1",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Notch1",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Cd34",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Gata1",  text_by = "label", text_colour = "black"),
  plotReducedDim(sce.filter.stem, dimred="UMAP", colour_by = "Gata2",  text_by = "label", text_colour = "black"),
  ncol=3)

#DE analysis
Sum_all <- aggregateAcrossCells(sce.filter.stem, id=colData(sce.filter.stem)[,c("Celltype","Genotype")])
Sum_CLP <- Sum_all[,which(Sum_all$Celltype =="CLP")]
DE_CLP <- DGEList(counts(Sum_CLP), samples = colData(Sum_CLP))
DE_CLP <- DE_CLP[Sum_CLP$ncells > 10,]
keep <- filterByExpr(DE_CLP, group= Sum_CLP$Genotype)
DE_CLP <- DE_CLP[keep,]
DE_CLP <- calcNormFactors(DE_CLP)
par(mfrow=c(2,2))
for (i in seq_len(ncol(DE_CLP))) {plotMD(DE_CLP, column=i)}
plotMDS(cpm(DE_CLP, log=TRUE), col=ifelse(DE_CLP$samples$condition, "red", "blue"))
DE_CLP <- estimateDisp(DE_CLP)
summary(DE_CLP$trended.dispersion)
design <- model.matrix(~factor(condition), DE_CLP$samples)
DE_CLP <- estimateDisp(DE_CLP, design)
summary(DE_CLP$trended.dispersion)
plotBCV(DE_CLP)
fit <- glmQLFit(DE_CLP, design, robust= TRUE)
plotQLDisp(fit)
res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))
nrow(res$table)
Pval_CLP <- as.data.frame(aaa$table)
Pval_CLP$log2FoldChange <-as.numeric(Pval_CLP$logFC)*(-1)
Pval_CLP <- merge(Pval_CLP, Ensembl.id.list, by=0)
Pval_CLP <- Pval_CLP[order(Pval_CLP$FDR),]
write.csv(Pval_CLP, "~/Desktop/scRNA-seq/Pval_CLP.csv")
EnhancedVolcano(Pval_CLP, lab = NA, x='log2FoldChange',y='FDR',
                title = 'CLP: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','darkblue'), colAlpha = 0.5,
                FCcutoff = 1.2, legendPosition = 'right') + xlim(-7.5,7.5) + ylim(0,12)
Genelist_CLP <- Pval_CLP$log2FoldChange
names(Genelist_CLP) <- Pval_CLP$Gene.stable.ID
Genelist_CLP <- na.omit(Genelist_CLP)
Genelist_CLP <- sort(Genelist_CLP, decreasing = TRUE)
Gse_CLP <- gseGO(geneList = Genelist_CLP, ont = "BP", keyType = "ENSEMBL", 
                        minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = TRUE,
                        OrgDb = "org.Mm.eg.db",pAdjustMethod = "none")
Gsea_CLP <- as.data.frame(Gse_CLP)
write.csv(Gsea_CLP, "~/Desktop/scRNA-seq/Gsea_CLP.csv")
require(DOSE)
dotplot(Gse_CLP, font.size=8, showCategory=8, split=".sign") + facet_grid(.~.sign) + ggtitle("CLP: Mut vs WT")
gseaplot2(Gse_CLP, geneSetID = 159, title = Gse_CLP$Description[159])
Gse_CLP$NES[159]
Gse_CLP$p.adjust[159]
#
Sum_CMP <- Sum_all[,which(Sum_all$Celltype =="CMP")]
DE_CMP <- DGEList(counts(Sum_CMP), samples = colData(Sum_CMP))
DE_CMP <- DE_CMP[Sum_CMP$ncells > 10,]
keep <- filterByExpr(DE_CMP, group= Sum_CMP$Genotype)
DE_CMP <- DE_CMP[keep,]
DE_CMP <- calcNormFactors(DE_CMP)
par(mfrow=c(2,2))
for (i in seq_len(ncol(DE_CMP))) {plotMD(DE_CMP, column=i)}
plotMDS(cpm(DE_CMP, log=TRUE), col=ifelse(DE_CMP$samples$condition, "red", "blue"))
DE_CMP <- estimateDisp(DE_CMP)
summary(DE_CMP$trended.dispersion)
design <- model.matrix(~factor(condition), DE_CMP$samples)
DE_CMP <- estimateDisp(DE_CMP, design)
summary(DE_CMP$trended.dispersion)
plotBCV(DE_CMP)
fit <- glmQLFit(DE_CMP, design, robust= TRUE)
plotQLDisp(fit)
res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))
nrow(res$table)
aaa<- topTags(res, n=12740)
Pval_CMP <- as.data.frame(aaa$table)
Pval_CMP$log2FoldChange <-as.numeric(Pval_CMP$logFC)*(-1)

Genelist_CMP <- Pval_CMP$log2FoldChange
names(Genelist_CMP) <- Pval_CMP$Gene.stable.ID
Genelist_CMP <- na.omit(Genelist_CMP)
Genelist_CMP <- sort(Genelist_CMP, decreasing = TRUE)
Gse_CMP <- gseGO(geneList = Genelist_CMP, ont = "BP", keyType = "ENSEMBL", 
                 minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.9, verbose = TRUE,
                 OrgDb = "org.Mm.eg.db",pAdjustMethod = "none")
Gsea_CMP <- as.data.frame(Gse_CMP)
write.csv(Gsea_CMP, "~/Desktop/scRNA-seq/Gsea_CMP.csv")
require(DOSE)
dotplot(Gse_CMP, font.size=8, showCategory=8, split=".sign") + facet_grid(.~.sign) + ggtitle("CLP: Mut vs WT")
Gse_CMP$NES[716]
Gse_CMP$p.adjust[716]

#Differential abundance
abundance <- table(sce.filter.stem$Celltype, sce.filter.stem$Genotype)
abundance <- unclass(abundance)
extra.info <- colData(sce.filter.stem)[match(colnames(abundance), sce.filter.stem$Genotype),]
extra.info$pred.immgen.fine<- NULL
extra.info$pred.immgen.main<- NULL
DE_ab <- DGEList(abundance, samples= extra.info)
keep <- filterByExpr(DE_ab, group= DE_ab$samples$condition)
DE_ab <- DE_ab[keep,]
design <- model.matrix(~factor(condition), DE_ab$samples)
DE_ab <- estimateDisp(DE_ab, design)
plotBCV(DE_ab)
fit.ab <- glmQLFit(DE_ab, design, robust=TRUE, abundance.trend=FALSE)
plotQLDisp(fit.ab, cex=1)
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)
write.csv(topTags(res), "~/Desktop/scRNA-seq/scRNAseq_DA.csv")

#Convert to seurat object
seurat.filter.stem <- sce.filter.stem
colnames(seurat.filter.stem)<- c(1:37685)
rownames(seurat.filter.stem)<- aaa@rowRanges$gene_id
reducedDim(seurat.filter.stem,"PCA.elbow") <- NULL
reducedDim(seurat.filter.stem,"PCA") <- NULL
seurat.filter.stem$pred.immgen.fine<-NULL
seurat.filter.stem$pred.immgen.main<-NULL
sce.seurat <- as.Seurat(seurat.filter.stem)
DimPlot(sce.seurat, reduction = "UMAP", group.by = "Celltype", label = FALSE ,pt.size= 1,
        cols = c("deeppink2","chartreuse1","blue3","darkturquoise","darkgreen","red","blueviolet","darkgrey","orange"))
FeaturePlot(sce.seurat, features = c("ENSMUSG00000005672","ENSMUSG00000075602","ENSMUSG00000016494","ENSMUSG00000015316" ) ,ncol=2)
FeaturePlot(sce.seurat, features = c("ENSMUSG00000015355","ENSMUSG00000002111","ENSMUSG00000034664","ENSMUSG00000042817" ) ,ncol=2)
FeaturePlot(sce.seurat, features = c("ENSMUSG00000054626"))
FeaturePlot(object, features = "ENSMUSG00000117465", split.by="condition", keep.scale = "all")
d2 <- as.data.frame(c("ENSMUSG00000117465","ENSMUSG00000024107","ENSMUSG00000028023","ENSMUSG00000022037","ENSMUSG00000032125"))
colnames(d2)<-"Ensemble.id"
object<- AddModuleScore(object = sce.seurat, features = d2, name = "lymphoid_score")
FeaturePlot(object, features = "lymphoid_score1", split.by="condition", keep.scale = "all", ncol = 1, max.cutoff = 0.4, min.cutoff = 0.000000001)
