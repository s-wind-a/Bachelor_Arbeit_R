library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("EnhancedVolcano")
library("rtracklayer")
library("ggpointdensity")
library("viridis")
library('biomaRt')
library("gridExtra")
library("grid")
library("UpSetR")
library("ComplexHeatmap")
library("egg")
library("matrixStats")
library("clusterProfiler")
library("enrichplot")
library("biomaRt")
library("rtracklayer")
library("org.Hs.eg.db")
library("ReactomePA")
library("eulerr")
library("GenomicFeatures")
library("rtracklayer")
library("GDCRNATools")
library("TCGAbiolinks")
library("SummarizedExperiment")
library("ComplexHeatmap")
library("UpSetR")

###########################################################################

###################### Functions

########### Functions for DESeq2 analysis

create_df_for_plot <- function(x){
  x$regulation <- "keine"
  x[which(x$log2FoldChange < -1 & x$padj < 0.05),]$regulation <- "2x runter"
  x[which(x$log2FoldChange > 1 & x$padj < 0.05),]$regulation <- "2x hoch"
  x$Expression <- log2(x$baseMean)
  return(x)
}

create_df_for_labels <- function(x){
  tab_ig <- x["ENSG00000124193",]
  tab_top_ig <- x[which(x$padj<0.05),]
  tab_top_ig <- tab_top_ig[order(tab_top_ig$padj),]
  tab_top_ig_up <- head(tab_top_ig[which(tab_top_ig$log2FoldChange>0),], 10)
  tab_top_ig_down <- head(tab_top_ig[which(tab_top_ig$log2FoldChange<0),], 10)
  tab_top_ig <- as.data.frame(rbind(tab_top_ig_down, tab_top_ig_up, tab_ig))
  symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(rownames(tab_top_ig), 15),mart= mart)
  tab_top_ig <- tab_top_ig[order(rownames(tab_top_ig)),]
  tab_top_ig$symbol <- symbols$hgnc_symbol
  return(tab_top_ig)
}

my_volcano_plot <- function(results, ig){
  ggplot(results[which(results$Expression>0),]) +
    geom_point(aes(y=-log10(padj), x=log2FoldChange, color = regulation)) +
    scale_color_manual(values = c("2x runter" = "midnightblue", "keine" = "lightyellow2", "2x hoch" = "violetred4")) +
    geom_vline(xintercept = -1:1, color = "grey44") +
    #    geom_vline(xintercept = 0.585, color = "grey44") +
    #    geom_vline(xintercept = -0.585, color = "grey44") +
    geom_hline(yintercept = 1.3, color = "grey44") +
    theme_minimal() + geom_label_repel(ig, 
                                       mapping = aes(label = symbol,
                                                     y = -log10(padj),
                                                     x = log2FoldChange),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 10, 
                                       fontface = "bold",
                                       size = 10/3)
}

my_ma_plot <- function(results, ig){
  ggplot(results[which(results$Expression>0),]) +
    ggtitle("Effekt der Hypoxie") +
    geom_point(aes(x=Expression, y=log2FoldChange, color = regulation)) +
    scale_color_manual(values = c("2x runter" = "midnightblue", "keine" = "lightyellow2", "2x hoch" = "violetred4")) +
    geom_hline(yintercept = -1:1, color = "grey44") +
    geom_hline(yintercept = 0.585, color = "grey44") +
    geom_hline(yintercept = -0.585, color = "grey44") +
    theme_minimal() + geom_label_repel(ig, 
                                       mapping = aes(label = symbol,
                                                     x = log2(baseMean),
                                                     y = log2FoldChange),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 1, 
                                       fontface = "bold",
                                       size = 10/3)
}

write_most_significant_genes <- function(reslfc, name){
  resOrdered <- reslfc[order(reslfc$padj),]
  resSig <- subset(resOrdered, padj < 0.05)
  symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(rownames(resSig), 15),mart= mart)
  rownames(symbols) <- paste(symbols$ensembl_gene_id)
  symbols <- symbols[strtrim(rownames(resSig), 15),]
  write.csv(symbols, paste0("D:/symbols_",name,".csv"))
  write.csv(resSig, paste0("D:/sig_",name,".csv"))
  return(resSig)
}

########### Function for visualisation of clusterProfiler analysis

my_dotplot <- function(x) {
  dotplot(x, showCategory=3) + 
    scale_color_viridis(option = "plasma", begin=0, end=0.75, direction = -1, guide=guide_colorbar(reverse=TRUE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(~ othergroup, scales = "free_x") +
    theme(strip.text.x = element_blank(),legend.text=element_text(size=10),
          legend.title=element_text(size=10))
}

########### Functions for analysis of MAJIQ results

create_igs_table_for_scatter <- function(scatter_table, gene_list){
  igs_table <- scatter_table[gene_list,]
  symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(rownames(igs_table), 15),mart= mart)
  igs_table$ensembl <- strtrim(rownames(igs_table),15)
  igs_table <- igs_table[!duplicated(igs_table$ensembl),]
  igs_table$true <- igs_table$ensembl%in%symbols$ensembl_gene_id
  igs_table <- igs_table[!duplicated(igs_table$ensembl),]
  igs_table <- igs_table[order(rownames(igs_table)),]
  igs_table$ensembl_2 <- symbols$ensembl_gene_id
  igs_table$symbol <- symbols$hgnc_symbol
  return(igs_table)
}

create_df_for_labels_splicing <- function(result, genes){
  tab_top_ig <- result[genes$V1,]
  symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(rownames(tab_top_ig), 15),mart= mart)
  tab_top_ig <- tab_top_ig[order(rownames(tab_top_ig)),]
  tab_top_ig$symbol <- symbols$hgnc_symbol
  return(tab_top_ig)
}

`%notin%` <- Negate(`%in%`)

########### Functions for comparison of results and data from TCGA

make_tcga_tpms <- function(x){
  
  CancerProject <- paste0("TCGA-", x)
  DataDirectory <- paste0("C:/Users/Ana/Desktop/3_counts/KICH", x)
  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")
  
  query <- GDCquery(project = CancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  samplesDown <- getResults(query,cols=c("cases"))
  
  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                    typesample = "TP")
  dataSmTP_short <- dataSmTP

  queryDown <- GDCquery(project = CancerProject, 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", 
                        barcode = c(dataSmTP_short))
  
  GDCdownload(query = queryDown,
              directory = DataDirectory)
  
  dataPrep <- GDCprepare(query = queryDown, 
                         save = TRUE, 
                         directory =  DataDirectory,
                         save.filename = FileNameData)
  
  dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep,
                                        #  cor.cut = 0.6,
                                        datatype = "HTSeq - Counts") 
  return(dataPrep)
}

count_tpm_2 <- function(count, gene_sizes){
  tpm_table <- as.data.frame(cbind(rownames(count), count[,n]))
  colnames(tpm_table) <- c("V1", "V2")
  tpm_table <- tpm_table[which(tpm_table$V1%in%strtrim(rownames(gene_sizes), 15)),]
  exons <- gene_sizes[tpm_table$V1,]
  tpm_table <- cbind(tpm_table, exons)
  tpm_table$exons_kb <- tpm_table$exons/1000
  tpm_table$rpk <- as.numeric(tpm_table$V2)/as.numeric(tpm_table$exons_kb)
  tpm_table$pm <- sum(tpm_table$rpk, na.rm=TRUE)/1000000
  tpm_table$tpm <- tpm_table$rpk/tpm_table$pm
  return(tpm_table$tpm)
}

make_table_for_scatters_tcga <- function(tpms) {
  correlation_srsf6 <- as.data.frame(cbind(SRSF6=c(log2(tpms["ENSG00000124193",]+1)),
                                           SRSF4=c(log2(tpms["ENSG00000116350",]+1)),
                                           MALAT1=c(log2(tpms["ENSG00000251562",]+1)),
                                           SRSF3=c(log2(tpms["ENSG00000112081",]+1)),
                                           SRSF7=c(log2(tpms["ENSG00000115875",]+1)),
                                           CHAF1A=c(log2(tpms["ENSG00000198625",]+1))))
  return(correlation_srsf6)
}

plot_green_tcga <- function(table, xx, yy){ # before I decided to do this plot blue it was green
  scatter <- ggplot(table, aes(x=unlist(xx), y=unlist(yy))) +
    geom_pointdensity() +
    geom_smooth(method = "lm", se=FALSE, colour="darkblue") +
    scale_color_viridis(option = "plasma", "neighbors") +
    theme_minimal() +
    scale_y_continuous(guide = guide_axis(position = "right"), limits = c(2, 13)) +
    xlim(c(4.5, 9)) +
    theme(panel.background = element_rect(fill = "#E5F7FF",
                                          colour = "#E5F7FF",
                                          size = 0.5, linetype = "solid")) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank()) +
    theme(legend.position = "none")
  return(scatter)
}

plot_red_tcga <- function(table, xx, yy){ # before I decided to do this plot yellow it was red
  scatter <- ggplot(table, aes(x=unlist(xx), y=unlist(yy))) +
    geom_pointdensity() +
    scale_color_viridis(option = "plasma", "neighbors") +
    geom_smooth(method = "lm", se=FALSE, colour="darkblue") +
    theme_minimal() +
    scale_y_continuous(guide = guide_axis(position = "right"), limits = c(2, 13)) +
    xlim(c(4.5, 9)) +
    theme(panel.background = element_rect(fill = "#FFFFDF",
                                          colour = "#FFFFDF",
                                          size = 0.5, linetype = "solid")) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
    theme(axis.title.y=element_blank()) +
    theme(legend.position = "none")
  return(scatter)
}

########### Function for findeng TPM values from HTSeq count files

count_tpm <- function(count, gene_sizes){
  tpm_table <- count
  tpm_table <- tpm_table[which(tpm_table$V1%in%rownames(gene_sizes)),]
  exons <- gene_sizes[tpm_table$V1,]
  tpm_table <- cbind(tpm_table, exons)
  tpm_table$exons_kb <- tpm_table$exons/1000
  tpm_table$rpk <- tpm_table$V2/tpm_table$exons_kb
  tpm_table$pm <- sum(tpm_table$rpk)/1000000
  tpm_table$tpm <- tpm_table$rpk/tpm_table$pm
  return(tpm_table)
}

###########################################################################

###################### Differential gene expression analysis

########### Loading the htseq-count files to the environment and preparing them

directory_all <- "D:/3_counts/all"
directory_oe <- "D:/3_counts/oe"

sampleFiles_all <- grep("file",list.files(directory_all),value=TRUE)
sampleFiles_oe <- grep("file",list.files(directory_oe),value=TRUE)

sampleCondition_all <- c('ne_hyp', 'ne_hyp','ne_norm','ne_norm','ne_norm','oe_hyp','oe_hyp','oe_norm', 'oe_norm')
sampleCondition_oe <- c('oe_hyp','oe_hyp','oe_norm', 'oe_norm')

sampleTable_all <- data.frame(sampleName = sampleFiles_all,
                              fileName = sampleFiles_all,
                              condition = sampleCondition_all)
sampleTable_oe <- data.frame(sampleName = sampleFiles_oe,
                             fileName = sampleFiles_oe,
                             condition = sampleCondition_oe)

ddsHTSeq_all <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_all,
                                           directory = directory_all,
                                           design= ~ condition)
ddsHTSeq_oe <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_oe,
                                          directory = directory_oe,
                                          design= ~ condition)

ddsHTSeq_all$condition <- factor(ddsHTSeq_all$condition, levels = c("ne_hyp","ne_norm", "oe_hyp", "oe_norm"))
ddsHTSeq_oe$condition <- factor(ddsHTSeq_oe$condition, levels = c("oe_hyp", "oe_norm"))

ddsHTSeq_all$condition <- relevel(ddsHTSeq_all$condition, ref = "ne_norm")
ddsHTSeq_oe$condition <- relevel(ddsHTSeq_oe$condition, ref = "oe_norm")

ddsHTSeq_all <- DESeq(ddsHTSeq_all)
ddsHTSeq_oe <- DESeq(ddsHTSeq_oe)

########### Principal component analysis

vsd_all <- vst(ddsHTSeq_all, blind=FALSE)
head(assay(vsd_all), 3)
pca_all <- plotPCA(vsd_all) +
  scale_color_viridis(discrete = TRUE, option = "viridis")+
  xlab("pc1 varianz")+
  ylab("pc2 varianz")+
  theme_minimal()

rv <- rowVars(assay(vsd_all))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd_all)[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:9)
colnames(scree_plot)<-c("varianz","komponentnummer")

scree <- ggplot(scree_plot, mapping=aes(x=komponentnummer, y=varianz*100))+
  ylab("varianz %")+
  geom_bar(stat="identity", fill = "darkblue")+
  scale_x_continuous(breaks = seq(0, 9, by = 1))+
  theme_minimal()

sampleDists <- dist(t(assay(vsd_all)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="-")
colnames(sampleDistMatrix) <- NULL
hm_all <- pheatmap(sampleDistMatrix,
                   clustering_distance_rows=sampleDists,
                   clustering_distance_cols=sampleDists,
                   color = inferno(10))

########### Differential gene expression analysis

res_all_ne <- results(ddsHTSeq_all, contrast = c("condition", "ne_hyp", "ne_norm"))
res_all_oe_vs_ne <- results(ddsHTSeq_all, contrast = c("condition", "oe_norm", "ne_norm"))
res_all_ne_norm_vs_oe_hyp <- results(ddsHTSeq_all, contrast = c("condition", "oe_hyp", "ne_norm"))

res_oe <- results(ddsHTSeq_oe, contrast = c("condition", "oe_hyp", "oe_norm"))

res_table_ne <- as.data.frame(res_all_ne)
res_table_oe_vs_ne <- as.data.frame(res_all_oe_vs_ne)
res_table_ne_norm_vs_oe_hyp <- as.data.frame(res_all_ne_norm_vs_oe_hyp)

res_table_oe <- as.data.frame(res_oe)

res_table_ne <- create_df_for_plot(res_table_ne)
res_table_oe_vs_ne <- create_df_for_plot(res_table_oe_vs_ne)
res_table_ne_norm_vs_oe_hyp <- create_df_for_plot(res_table_ne_norm_vs_oe_hyp)
res_table_oe <- create_df_for_plot(res_table_oe)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

top_ne_ig <- create_df_for_labels(res_table_ne)
top_oe_vs_ne_ig <- create_df_for_labels(res_table_oe_vs_ne)
top_ne_norm_vs_oe_hyp_ig <- create_df_for_labels(res_table_ne_norm_vs_oe_hyp)
top_oe_ig <- create_df_for_labels(res_table_oe)

v1 <- my_volcano_plot(res_table_ne, top_ne_ig) + theme(legend.position = "none") + xlim(c(-25,25))
v2 <- my_volcano_plot(res_table_oe_vs_ne, top_oe_vs_ne_ig) + theme(legend.position = "none") + xlim(c(-25,25))
v3 <- my_volcano_plot(res_table_ne_norm_vs_oe_hyp, top_ne_norm_vs_oe_hyp_ig) + theme(legend.position = "none") + xlim(c(-25,25))
v4 <- my_volcano_plot(res_table_oe, top_oe_ig) + xlim(c(-25,25))
ggarrange(v1, v2, v3, v4, ncol = 2, labels=LETTERS[1:4], label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))

resLFC_all_ne <- lfcShrink(ddsHTSeq_all, coef="condition_ne_hyp_vs_ne_norm", type="apeglm")
resLFC_all_oe_vs_ne <- lfcShrink(ddsHTSeq_all, coef="condition_oe_norm_vs_ne_norm", type="apeglm")
resLFC_all_ne_norm_vs_oe_hyp <- lfcShrink(ddsHTSeq_all, coef="condition_oe_hyp_vs_ne_norm", type="apeglm")
resLFC_oe <- lfcShrink(ddsHTSeq_oe, coef="condition_oe_hyp_vs_oe_norm", type="apeglm")

resLFC_table_ne <- as.data.frame(resLFC_all_ne)
resLFC_table_oe_vs_ne <- as.data.frame(resLFC_all_oe_vs_ne)
resLFC_table_ne_norm_vs_oe_hyp <- as.data.frame(resLFC_all_ne_norm_vs_oe_hyp)
resLFC_table_oe <- as.data.frame(resLFC_oe)

resSig_ne <- write_most_significant_genes(resLFC_table_ne, "ne")
resSig_oe_vs_ne <- write_most_significant_genes(resLFC_table_oe_vs_ne, "oe_vs_ne")
resSig_ne_norm_vs_oe_hyp <- write_most_significant_genes(resLFC_table_ne_norm_vs_oe_hyp, "ne_norm_vs_oe_hyp")
resSig_oe <- write_most_significant_genes(resLFC_table_oe, "oe")

resLFC_table_ne <- create_df_for_plot(resLFC_table_ne)
resLFC_table_oe_vs_ne <- create_df_for_plot(resLFC_table_oe_vs_ne)
resLFC_table_ne_norm_vs_oe_hyp <- create_df_for_plot(resLFC_table_ne_norm_vs_oe_hyp)
resLFC_table_oe <- create_df_for_plot(resLFC_table_oe)

tab_top_ne_ig <- create_df_for_labels(resLFC_table_ne)
tab_top_oe_vs_ne_ig <- create_df_for_labels(resLFC_table_oe_vs_ne)
top_ne_norm_vs_oe_hyp_ig <- create_df_for_labels(resLFC_table_ne_norm_vs_oe_hyp)
top_oe_ig <- create_df_for_labels(resLFC_table_oe)

ma1 <- my_ma_plot(resLFC_table_ne, tab_top_ne_ig)
ma2 <- my_ma_plot(resLFC_table_oe_vs_ne, tab_top_oe_vs_ne_ig)
ma3 <- my_ma_plot(resLFC_table_ne_norm_vs_oe_hyp, top_ne_norm_vs_oe_hyp_ig)
ma4 <- my_ma_plot(resLFC_table_oe, top_oe_ig)

scatter_table_oe_vs_ne <- cbind(resLFC_table_oe, resLFC_table_ne)
colnames(scatter_table_oe_vs_ne) <- c("baseMean_oe","log2FoldChange_oe","lfcSE_oe","pvalue_oe","padj_oe",
                                      "regulation_oe","Expression_oe","baseMean_ne","log2FoldChange_ne","lfcSE_ne",
                                      "pvalue_ne","padj_ne","regulation_ne","Expression_ne")
tab_top_ne_ig <- tab_top_ne_ig[rownames(tab_top_ne_ig)!="ENSG00000124193.16",]
top_oe_ig <- top_oe_ig[rownames(top_oe_ig)!="ENSG00000124193.16",]
scatter_table_oe_vs_ne_ig <- rbind(scatter_table_oe_vs_ne["ENSG00000124193.16",],
                                   scatter_table_oe_vs_ne[rownames(tab_top_ne_ig),],
                                   scatter_table_oe_vs_ne[rownames(top_oe_ig),])
symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(rownames(scatter_table_oe_vs_ne_ig), 15),mart= mart)
scatter_table_oe_vs_ne_ig <- scatter_table_oe_vs_ne_ig[!duplicated(rownames(scatter_table_oe_vs_ne_ig)),]
scatter_table_oe_vs_ne_ig$ensembl <- strtrim(rownames(scatter_table_oe_vs_ne_ig),15)
scatter_table_oe_vs_ne_ig$true <- scatter_table_oe_vs_ne_ig$ensembl%in%symbols$ensembl_gene_id
scatter_table_oe_vs_ne_ig <- scatter_table_oe_vs_ne_ig[!duplicated(scatter_table_oe_vs_ne_ig$ensembl),]
scatter_table_oe_vs_ne_ig <- scatter_table_oe_vs_ne_ig[order(rownames(scatter_table_oe_vs_ne_ig)),]
scatter_table_oe_vs_ne_ig$ensembl_2 <- symbols$ensembl_gene_id
scatter_table_oe_vs_ne_ig$symbol <- symbols$hgnc_symbol

ggplot(scatter_table_oe_vs_ne, aes(x=log2FoldChange_ne, y=log2FoldChange_oe)) +
  geom_pointdensity() +
  scale_color_viridis(option = "plasma", "neighbors") +
  theme_minimal() + geom_label_repel(scatter_table_oe_vs_ne_ig, 
                                     mapping = aes(label = symbol,
                                                   y = log2FoldChange_oe,
                                                   x = log2FoldChange_ne),
                                     box.padding = unit(0.35, "lines"),
                                     point.padding = unit(0.3, "lines"),
                                     force = 1, 
                                     fontface = "bold",
                                     size = 10/3)

########### Analysis of trends of all gene expression levels

violin <- data.frame(
  log2FoldChange = c(resSig_ne$log2FoldChange, resSig_oe_vs_ne$log2FoldChange, resSig_ne_norm_vs_oe_hyp$log2FoldChange, resSig_oe$log2FoldChange),
  Kondition = c(rep("a", length(rownames(resSig_ne))), rep("b",length(rownames(resSig_oe_vs_ne))),
                rep("c", length(rownames(resSig_ne_norm_vs_oe_hyp))), rep("d", length(rownames(resSig_oe))))
)
vio1 <- ggplot(violin, aes(x=Kondition, y=log2FoldChange, fill=Kondition)) + 
  scale_fill_viridis(discrete=T, name="", option = "plasma") +
  geom_violin() +
  geom_boxplot(width=0.5, fill=NA, colour="black") +
  ylim(c(-20,20)) +
  theme_minimal()+ theme(legend.position = "none")

violin2 <- data.frame(
  log2FoldChange = c(res_all_ne$log2FoldChange, res_all_oe_vs_ne$log2FoldChange, res_all_ne_norm_vs_oe_hyp$log2FoldChange, res_oe$log2FoldChange),
  Kondition = c(rep("a", length(rownames(res_all_ne))), rep("b",length(rownames(res_all_oe_vs_ne))),
                rep("c", length(rownames(res_all_ne_norm_vs_oe_hyp))), rep("d", length(rownames(res_oe))))
)
vio2 <- ggplot(violin2, aes(x=Kondition, y=log2FoldChange, fill=Kondition)) + 
  scale_fill_viridis(discrete=T, name="", option = "plasma") +
  geom_violin() +
  geom_boxplot(width=0.5, fill=NA, colour="black") +
  ylim(c(-20,20)) +
  theme_minimal() + theme(legend.position = "none")

pcas <- cowplot::plot_grid(pca_all, scree, ncol=2, labels=LETTERS[1:2])
vios <- ggarrange(vio2, vio1, ncol=2, labels=LETTERS[3:4], label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
grid.arrange(pcas, vios, ncol = 1, heights = c(3,5))

########### Comparison of differentialy expressed genes

venn1 <- list(a=rownames(resSig_ne[which(resSig_ne$log2FoldChange > 0),]),
              b=rownames(resSig_oe_vs_ne[which(resSig_oe_vs_ne$log2FoldChange > 0),]),
              c=rownames(resSig_ne_norm_vs_oe_hyp[which(resSig_ne_norm_vs_oe_hyp$log2FoldChange > 0),]),
              d=rownames(resSig_oe[which(resSig_oe$log2FoldChange > 0),]))
venn2 <- list(a=rownames(resSig_ne[which(resSig_ne$log2FoldChange < 0),]),
              b=rownames(resSig_oe_vs_ne[which(resSig_oe_vs_ne$log2FoldChange < 0),]),
              c=rownames(resSig_ne_norm_vs_oe_hyp[which(resSig_ne_norm_vs_oe_hyp$log2FoldChange < 0),]),
              d=rownames(resSig_oe[which(resSig_oe$log2FoldChange < 0),]))

upset1 <- list_to_matrix(venn1)
upset2 <-list_to_matrix(venn2)
upset1_m = make_comb_mat(upset1)
upset2_m = make_comb_mat(upset2)

u_my1 <- UpSet(upset1_m, top_annotation = HeatmapAnnotation(
  "Überlappende\nGene" = anno_barplot(comb_size(upset1_m),
                                      border = FALSE, 
                                      gp = gpar(fill = plasma(4)[comb_degree(upset1_m)]),
                                      height = unit(6, "cm"), ylim = c(0, 3000)),
  annotation_name_side = "left", annotation_name_rot = 0),
  left_annotation = rowAnnotation(
    "Hochregulierte Gene" = anno_barplot(set_size(upset1_m), 
                                         axis_param = list(direction = "reverse"),
                                         border = FALSE, 
                                         gp = gpar(fill = plasma(1)), 
                                         width = unit(4, "cm"), ylim = c(0, 8000),
                                         set_size.show = TRUE)),
  right_annotation = NULL)

u_my2 <- UpSet(upset2_m, top_annotation = HeatmapAnnotation(
  "Überlappende\nGene" = anno_barplot(comb_size(upset2_m),
                                      border = FALSE, 
                                      gp = gpar(fill = plasma(4)[comb_degree(upset2_m)]),
                                      height = unit(6, "cm"), ylim = c(0, 3000)),
  annotation_name_side = "left", annotation_name_rot = 0),
  left_annotation = rowAnnotation(
    "Herunterregulierte Gene" = anno_barplot(set_size(upset2_m), 
                                             axis_param = list(direction = "reverse"),
                                             border = FALSE, 
                                             gp = gpar(fill = plasma(1)), 
                                             width = unit(4, "cm"), ylim = c(0, 8000))), right_annotation = NULL)
u_my1+u_my2

###########################################################################

###################### Gene enrichment analysis

########### Preparing for analysis

go_ne <- resSig_ne
go_ne$name <- strtrim(rownames(resSig_ne), 15)
go_ne$oreg <- "1"
go_ne_up <- go_ne[go_ne$log2FoldChange > 1,]
go_ne_up <- go_ne_up[order(go_ne_up$log2FoldChange, decreasing = TRUE),]
go_ne_down <- go_ne[go_ne$log2FoldChange < 1,]
go_ne_down <- go_ne_down[order(go_ne_down$log2FoldChange, decreasing = FALSE),]
go_oe <- resSig_oe_vs_ne
go_oe$name <- strtrim(rownames(resSig_oe_vs_ne), 15)
go_oe$oreg <- "2"
go_oe_up <- go_oe[go_oe$log2FoldChange > 1,]
go_oe_up <- go_oe_up[order(go_oe_up$log2FoldChange, decreasing = TRUE),]
go_oe_down <- go_oe[go_oe$log2FoldChange < 1,]
go_oe_down <- go_oe_down[order(go_oe_down$log2FoldChange, decreasing = FALSE),]
go_oe_hyp <- resSig_ne_norm_vs_oe_hyp
go_oe_hyp$name <- strtrim(rownames(resSig_ne_norm_vs_oe_hyp), 15)
go_oe_hyp$oreg <- "3"
go_oe_hyp_up <- go_oe_hyp[go_oe_hyp$log2FoldChange > 1,]
go_oe_hyp_up <- go_oe_hyp_up[order(go_oe_hyp_up$log2FoldChange, decreasing = TRUE),]
go_oe_hyp_down <- go_oe_hyp[go_oe_hyp$log2FoldChange < 1,]
go_oe_hyp_down <- go_oe_hyp_down[order(go_oe_hyp_down$log2FoldChange, decreasing = FALSE),]
go_oe_norm_vs_hyp <- resSig_oe
go_oe_norm_vs_hyp$name <- strtrim(rownames(resSig_oe), 15)
go_oe_norm_vs_hyp$oreg <- "4"
go_oe_norm_vs_hyp_up <- go_oe_norm_vs_hyp[go_oe_norm_vs_hyp$log2FoldChange > 1,]
go_oe_norm_vs_hyp_up <- go_oe_norm_vs_hyp_up[order(go_oe_norm_vs_hyp_up$log2FoldChange, decreasing = TRUE),]
go_oe_norm_vs_hyp_down <- go_oe_norm_vs_hyp[go_oe_norm_vs_hyp$log2FoldChange < 1,]
go_oe_norm_vs_hyp_down <- go_oe_norm_vs_hyp_down[order(go_oe_norm_vs_hyp_down$log2FoldChange, decreasing = FALSE),]
go_oe_hyp_up$reg <- go_oe_up$reg <- go_ne_up$reg <- go_oe_norm_vs_hyp_up$reg <- "a"
go_oe_hyp_down$reg <- go_ne_down$reg <- go_oe_down$reg <- go_oe_norm_vs_hyp_down$reg <- "b"

ez_ne_down <- bitr(geneID = go_ne_down$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_oe_down <- bitr(geneID = go_oe_down$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_oe_hyp_down <- bitr(geneID = go_oe_norm_vs_hyp_down$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_oe_norm_vs_hyp_down <- bitr(geneID = go_oe_norm_vs_hyp_down$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_ne_up <- bitr(geneID = go_ne_up$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_oe_up <- bitr(geneID = go_oe_up$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_oe_hyp_up <- bitr(geneID = go_oe_hyp_up$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_oe_norm_vs_hyp_up <- bitr(geneID = go_oe_norm_vs_hyp_up$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
ez_ne_down$reg <- ez_oe_down$reg <- ez_oe_hyp_down$reg <- ez_oe_norm_vs_hyp_down$reg <- "b"
ez_ne_up$reg <- ez_oe_up$reg <- ez_oe_hyp_up$reg <- ez_oe_norm_vs_hyp_up$reg <- "a"
ez_ne_down$oreg <- ez_ne_up$oreg <- "1"
ez_oe_down$oreg <- ez_oe_up$oreg <- "2"
ez_oe_hyp_down$oreg <- ez_oe_hyp_up$oreg <- "3"
ez_oe_norm_vs_hyp_down$oreg <- ez_oe_norm_vs_hyp_up$oreg <- "4"

godf_ez <- data.frame(Ensembl = c(ez_ne_up$ENTREZID, ez_ne_down$ENTREZID, ez_oe_up$ENTREZID, ez_oe_down$ENTREZID, ez_oe_hyp_up$ENTREZID,
                                  ez_oe_hyp_down$ENTREZID, ez_oe_norm_vs_hyp_up$ENTREZID, ez_oe_norm_vs_hyp_down$ENTREZID),
                      group = c(ez_ne_up$oreg, ez_ne_down$oreg, ez_oe_up$oreg, ez_oe_down$oreg, ez_oe_hyp_up$oreg,ez_oe_hyp_down$oreg,
                                ez_oe_norm_vs_hyp_up$oreg, ez_oe_norm_vs_hyp_down$oreg),
                      othergroup = c(ez_ne_up$reg, ez_ne_down$reg, ez_oe_up$reg, ez_oe_down$reg, ez_oe_hyp_up$reg,ez_oe_hyp_down$reg,
                                     ez_oe_norm_vs_hyp_up$reg, ez_oe_norm_vs_hyp_down$reg))

univ <- rbind(go_ne, go_oe, go_oe_hyp)
univ <- univ[order(univ$log2FoldChange, decreasing = TRUE),]
univ_ez <- bitr(geneID = univ$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)

# this universe used because I wanted to see only the regulated processes in one or more conditions to see more differences between them
# with a default universe less differences are detected and many processes are up- and downregulated in the same time in same kondition

########### Analysis and visualisation

compareBP_ez <- compareCluster(Ensembl~group+othergroup,
                               data          = godf_ez,
                               universe      = univ_ez$ENTREZID,
                               fun           = "enrichGO",
                               OrgDb         = org.Hs.eg.db,
                               keyType       = 'ENTREZID',
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)

bp <- my_dotplot(compareBP_ez)

compareCC_ez <- compareCluster(Ensembl~group+othergroup,
                               data          = godf_ez,
                               universe      = univ_ez$ENTREZID,
                               fun           = "enrichGO",
                               OrgDb         = org.Hs.eg.db,
                               keyType       = 'ENTREZID',
                               ont           = "CC",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)

cc <- my_dotplot(compareCC_ez)

compareMF_ez <- compareCluster(Ensembl~group+othergroup,
                               data          = godf_ez,
                               universe      = univ_ez$ENTREZID,
                               fun           = "enrichGO",
                               OrgDb         = org.Hs.eg.db,
                               keyType       = 'ENTREZID',
                               ont           = "MF",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)

mf <- my_dotplot(compareMF_ez)

compareKEGG <- compareCluster(Ensembl~group+othergroup,
                              data          = godf_ez,
                              universe      = univ_ez$ENTREZID,
                              fun           = "enrichKEGG",
                              organism     = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05)

kegg <- my_dotplot(compareKEGG)

comparePA <- compareCluster(Ensembl~group+othergroup,
                            data          = godf_ez,
                            universe      = univ_ez$ENTREZID,
                            fun           = "enrichPathway",
                            organism     = 'human',
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.05)

pa <- my_dotplot(comparePA)

ggarrange(bp,cc,mf, labels=LETTERS[1:3], label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
ggarrange(pa,kegg, labels=LETTERS[1:2], heights = c(5,4), label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))

###########################################################################

###################### Analysis of MAJIQ/Voila output

########### Prepare output for analysis and visualisation

deltapsi_ne_tsv <- read.table(file = 'D:/deltapsi_fin_strikt/wt_normoxia_wt_hypoxia.deltapsi.tsv', sep = '\t', header = TRUE)
deltapsi_oe_tsv <- read.table(file = 'D:/deltapsi_fin_strikt/oe_normoxia_oe_hypoxia.deltapsi.tsv', sep = '\t', header = TRUE)

deltapsi_ne <- read.table(file = 'D:/deltapsi/wt_normoxia_wt_hypoxia.txt', header = FALSE)
deltapsi_oe <- read.table(file = 'D:/deltapsi/oe_normoxia_oe_hypoxia.txt', header = FALSE)

genes_ne <- read.table(file = 'D:/deltapsi/genes_ne.txt', header = FALSE)
genes_oe <- read.table(file = 'D:/deltapsi/genes_oe.txt', header = FALSE)

########### Quantitative analysis of alternative spliced genes

eul <- list(a = deltapsi_ne$V1,
            b = deltapsi_oe$V1)

eul2 <- list(a = genes_ne$V1,
             b = genes_oe$V1)

eup1 <- plot(euler(eul, shape = "ellipse"), quantities = TRUE, fills = c("coral1","yellow", "#FFF3BE"))
eup2 <- plot(euler(eul2, shape = "ellipse"), quantities = TRUE, fills = c("coral1","yellow", "#FFF3BE"))

########### Quantitative analysis of alternative splicing events, interactive splice graphs used

five_prime_ne <- 25
five_prime_only_ne <- 3
five_prime_other_ne <- five_prime_ne - five_prime_only_ne
five_prime_ne_up <- 9
five_prime_ne_down <- 5
three_prime_ne <- 23
three_prime_only_ne <- 4
three_prime_other_ne <- three_prime_ne - three_prime_only_ne
three_prime_ne_up <- 6
three_prime_ne_down <- 3
exon_skipping_ne <- 124
exon_skipping_only_ne <- 44
exon_skipping_other_ne <- exon_skipping_ne - exon_skipping_only_ne
exon_skipping_ne_up <- 67
exon_skipping_ne_down <- 29
intron_retention_ne <- 34
intron_retention_only_ne <- 3
intron_retention_other_ne <- intron_retention_ne - intron_retention_only_ne
intron_retention_ne_up <- 9
intron_retention_ne_down <- 7
sum_ne <- five_prime_ne+three_prime_ne+exon_skipping_ne+intron_retention_ne

five_prime_oe <- 17
five_prime_only_oe <- 1
five_prime_other_oe <- five_prime_oe - five_prime_only_oe
five_prime_oe_down <- 4
five_prime_oe_up <- 5
three_prime_oe <- 9
three_prime_only_oe <- 1
three_prime_other_oe <- three_prime_oe - three_prime_only_oe
three_prime_oe_up <- 1
three_prime_oe_down <- 3
exon_skipping_oe <- 52
exon_skipping_only_oe <- 16
exon_skipping_other_oe <- exon_skipping_oe - exon_skipping_only_oe
exon_skipping_oe_up <- 26
exon_skipping_oe_down <- 11
intron_retention_oe <- 37
intron_retention_only_oe <- 12
intron_retention_other_oe <- intron_retention_oe - intron_retention_only_oe
intron_retention_oe_up <- 12
intron_retention_oe_down <- 8
sum_oe <- five_prime_oe+three_prime_oe+exon_skipping_oe+intron_retention_oe

splicing_reg_df <- data.frame(spleißereignis = c("5' ASS", "3' ASS","AE", "IR", "5' ASS", "3' ASS","AE", "IR",
                                                 "5' ASS", "3' ASS","AE", "IR", "5' ASS", "3' ASS","AE", "IR"),
                              regulation = c("hoch", "hoch", "hoch", "hoch", "runter", "runter", "runter", "runter",
                                             "hoch", "hoch", "hoch", "hoch", "runter", "runter", "runter", "runter"),
                              srsf6 = c("ne","ne","ne","ne","ne","ne","ne","ne",
                                        "oe","oe","oe","oe","oe","oe","oe","oe"),
                              zahl = c(five_prime_ne_up/sum_ne, three_prime_ne_up/sum_ne, exon_skipping_ne_up/sum_ne, intron_retention_ne_up/sum_ne,
                                       five_prime_ne_down/sum_ne, three_prime_ne_down/sum_ne, exon_skipping_ne_down/sum_ne, intron_retention_ne_down/sum_ne,
                                       five_prime_oe_up/sum_oe, three_prime_oe_up/sum_oe, exon_skipping_oe_up/sum_oe, intron_retention_oe_up/sum_oe,
                                       five_prime_oe_down/sum_oe, three_prime_oe_down/sum_oe, exon_skipping_oe_down/sum_oe, intron_retention_oe_down/sum_oe))

splicing_reg_df$regulation <- factor(splicing_reg_df$regulation, levels = c("hoch", "runter"))

ggplot(data=splicing_reg_df, aes(x=spleißereignis, y=zahl*100, fill=regulation)) +
  ylab("prozent")+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("gold", "dodgerblue4"))+
  facet_grid(~srsf6) +
  theme_bw()

out_df_ne <- data.frame(
  spleißereignis = c("5' ASS binär","5' ASS komplex", "3' ASS binär","3' ASS komplex", "AE binär",
                     "AE komplex","IR binär","IR komplex", "5' ASS", "3' ASS",
                     "AE", "IR"),
  value = c(as.numeric(five_prime_only_ne), as.numeric(five_prime_other_ne), as.numeric(three_prime_only_ne),as.numeric(three_prime_other_ne),
            as.numeric(exon_skipping_only_ne),as.numeric(exon_skipping_other_ne), as.numeric(intron_retention_only_ne),
            as.numeric(intron_retention_other_ne), as.numeric(five_prime_ne), as.numeric(three_prime_ne),
            as.numeric(exon_skipping_ne), as.numeric(intron_retention_ne)),
  float = c(rep(2,8), rep(1,4))
)

out_df_oe <- data.frame(
  spleißereignis = c("5' ASS binär","5' ASS komplex", "3' ASS binär","3' ASS komplex", "AE binär",
                     "AE komplex","IR binär","IR komplex", "5' ASS", "3' ASS",
                     "AE", "IR"),
  value = c(as.numeric(five_prime_only_oe), as.numeric(five_prime_other_oe), as.numeric(three_prime_only_oe),as.numeric(three_prime_other_oe),
            as.numeric(exon_skipping_only_oe),as.numeric(exon_skipping_other_oe), as.numeric(intron_retention_only_oe),
            as.numeric(intron_retention_other_oe), as.numeric(five_prime_oe), as.numeric(three_prime_oe),
            as.numeric(exon_skipping_oe), as.numeric(intron_retention_oe)),
  float = c(rep(2,8), rep(1,4))
)

pie1 <- ggplot(out_df_ne, aes(x = float,y = value,fill = spleißereignis)) + 
  geom_bar(width = 1, stat = "identity") + 
  scale_fill_manual(values = c("midnightblue", "#B2D5FF", "whitesmoke", "violetred4", "#9D75AF",
                               "#F8F3FF", "coral", "#67ABC6", "#FFEDE9", "yellow2", "#89C2AD", "#FCFEEC")) + 
  coord_polar("y", start = 2*pi) +
  theme_void() +
  theme(legend.position = 'none',legend.text=element_text(size=10),
        legend.title=element_text(size=10))


pie2 <- ggplot(out_df_oe, aes(x = float,y = value,fill = spleißereignis)) + 
  geom_bar(width = 1, stat = "identity") + 
  scale_fill_manual(values = c("midnightblue", "#B2D5FF", "whitesmoke", "violetred4", "#9D75AF",
                               "#F8F3FF", "coral", "#67ABC6", "#FFEDE9", "yellow2", "#89C2AD", "#FCFEEC"),
                    guide=guide_legend(nrow=2)) + 
  coord_polar("y", start = 2*pi) +
  theme_void() +
  theme(legend.position = "bottom",legend.text=element_text(size=10),
        legend.title=element_text(size=10))

legend <- cowplot::get_legend(pie2)

pie2 <- pie2+theme(legend.position = 'none')

pie <- cowplot::plot_grid(pie1,pie2, ncol=2, labels=LETTERS[1:2], rel_heights = c(3,5,1))
cowplot::plot_grid(pie, legend, ncol=1, rel_heights = c(9,1))
cowplot::plot_grid(eup1, eup2, labels=LETTERS[1:2], ncol=2)

########### Correlation of alternative splicing events and differential gene expression

symbols_ne <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(genes_ne$V1, 15),mart= mart)
symbols_oe <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(genes_oe$V1, 15),mart= mart)

write_xlsx(as.data.frame(symbols_ne), 'D:/deltapsi/symbols_ne.xlsx')
write_xlsx(as.data.frame(symbols_oe), 'D:/deltapsi/symbols_oe.xlsx')

genes_only_ne <- genes_ne[which(genes_ne$V1 %notin% genes_oe$V1),]
genes_only_oe <- genes_oe[which(genes_oe$V1 %notin% genes_ne$V1),]

eul_splicing_vs_deg <- list(a=resLFC_all_ne[genes_ne$V1,]$log2FoldChange,b=resLFC_oe[genes_oe$V1,]$log2FoldChange)

violin_splicing_vs_deg <- data.frame(
  log2FoldChange = c(resLFC_all_ne[genes_ne$V1,]$log2FoldChange, resLFC_oe[genes_oe$V1,]$log2FoldChange),
  Kondition = c(rep("a", length(rownames(resLFC_all_ne[genes_ne$V1,]))), rep("b", length(rownames(resLFC_oe[genes_oe$V1,])))))

violin_splicing <- ggplot(violin_splicing_vs_deg, aes(x=Kondition, y=log2FoldChange, fill=Kondition)) + 
  # ggtitle("C") +
  scale_fill_manual(values=c("coral1", "yellow")) +
  geom_violin() +
  geom_boxplot(width=0.5, fill=NA, colour="black") +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_minimal() +
  theme(legend.position = "none")

labels_ne <- create_df_for_labels_splicing(res_table_ne, genes_ne)
labels_oe <- create_df_for_labels_splicing(res_table_oe, genes_oe)
volcano_oe_splicing <- my_volcano_plot(res_table_oe, labels_oe[which(labels_oe$log2FoldChange > 1 | labels_oe$log2FoldChange < -1),]) + xlim(c(-10,25))
volcano_ne_splicing <- my_volcano_plot(res_table_ne, labels_ne[which(labels_ne$log2FoldChange > 1 | labels_ne$log2FoldChange < -1),]) + xlim(c(-10,25)) + theme(legend.position = "none")

genes_both <- genes_ne[which(genes_ne$V1 %in% genes_oe$V1),]

scatter_table_both_ig <- create_igs_table_for_scatter(scatter_table_oe_vs_ne, genes_both)

scatter_splicing_both <- ggplot(scatter_table_oe_vs_ne, aes(x=log2FoldChange_ne, y=log2FoldChange_oe)) +
  geom_pointdensity() +
  scale_color_viridis(option = "plasma", "neighbors") +
  theme_minimal() + geom_label_repel(scatter_table_both_ig, 
                                     mapping = aes(label = symbol,
                                                   y = log2FoldChange_oe,
                                                   x = log2FoldChange_ne),
                                     box.padding = unit(0.35, "lines"),
                                     point.padding = unit(0.3, "lines"),
                                     force = 1, 
                                     fontface = "bold",
                                     size = 10/3)

p1 <- ggarrange(volcano_ne_splicing, volcano_oe_splicing, labels=LETTERS[3:4], ncol=2, label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
p2 <- cowplot::plot_grid(violin_splicing, scatter_splicing_both, labels=LETTERS[5:6], ncol=2, rel_widths = c(1,3))
cowplot::plot_grid(p1,p2, ncol=1, rel_heights = c(5,4))

###########################################################################

###################### Comparison with TCGA

########### Prepare TPM-Values

dlbc_tpms_1 <- make_tcga_tpms("DLBC") 
dlbc_tpms <- as.data.frame(cbind(rownames(dlbc_tpms_1), rownames(dlbc_tpms_1)))
rownames(dlbc_tpms) <- rownames(dlbc_tpms_1)
for (i in 1:length(colnames(dlbc_tpms_1))) {dlbc_tpms[,i] <- NA}
colnames(dlbc_tpms) <- colnames(dlbc_tpms_1)
dlbc_tpms_1 <- as.data.frame(dlbc_tpms_1)
for (n in 1:length(colnames(dlbc_tpms_1))) {
  dlbc_tpms[,n] <- count_tpm_2(dlbc_tpms_1, exonic.gene.sizes)
  print(n)
}
gbm_tpms_1 <- make_tcga_tpms("GBM") 
gbm_tpms <- as.data.frame(cbind(rownames(gbm_tpms_1), rownames(gbm_tpms_1)))
rownames(gbm_tpms) <- rownames(gbm_tpms_1)
for (i in 1:length(colnames(gbm_tpms_1))) {gbm_tpms[,i] <- NA}
colnames(gbm_tpms) <- colnames(gbm_tpms_1)
gbm_tpms_1 <- as.data.frame(gbm_tpms_1)
for (n in 1:length(colnames(gbm_tpms_1))) {
  gbm_tpms[,n] <- count_tpm_2(gbm_tpms_1, exonic.gene.sizes)
  print(n)
}
lgg_tpms_1 <- make_tcga_tpms("LGG") 
lgg_tpms <- as.data.frame(cbind(rownames(lgg_tpms_1), rownames(lgg_tpms_1)))
rownames(lgg_tpms) <- rownames(lgg_tpms_1)
for (i in 1:length(colnames(lgg_tpms_1))) {lgg_tpms[,i] <- NA}
colnames(lgg_tpms) <- colnames(lgg_tpms_1)
lgg_tpms_1 <- as.data.frame(lgg_tpms_1)
for (n in 1:length(colnames(lgg_tpms_1))) {
  lgg_tpms[,n] <- count_tpm_2(lgg_tpms_1, exonic.gene.sizes)
  print(n)
}
thym_tpms_1 <- make_tcga_tpms("THYM") 
thym_tpms <- as.data.frame(cbind(rownames(thym_tpms_1), rownames(thym_tpms_1)))
rownames(thym_tpms) <- rownames(thym_tpms_1)
for (i in 1:length(colnames(thym_tpms_1))) {thym_tpms[,i] <- NA}
colnames(thym_tpms) <- colnames(thym_tpms_1)
thym_tpms_1 <- as.data.frame(thym_tpms_1)
for (n in 1:length(colnames(thym_tpms_1))) {
  thym_tpms[,n] <- count_tpm_2(thym_tpms_1, exonic.gene.sizes)
  print(n)
}
kich_tpms_1 <- make_tcga_tpms("KICH") 
kich_tpms <- as.data.frame(cbind(rownames(kich_tpms_1), rownames(kich_tpms_1)))
rownames(kich_tpms) <- rownames(kich_tpms_1)
for (i in 1:length(colnames(kich_tpms_1))) {kich_tpms[,i] <- NA}
colnames(kich_tpms) <- colnames(kich_tpms_1)
kich_tpms_1 <- as.data.frame(kich_tpms_1)
for (n in 1:length(colnames(kich_tpms_1))) {
  kich_tpms[,n] <- count_tpm_2(kich_tpms_1, exonic.gene.sizes)
  print(n)
}
ov_tpms_1 <- make_tcga_tpms("OV") 
ov_tpms <- as.data.frame(cbind(rownames(ov_tpms_1), rownames(ov_tpms_1)))
rownames(ov_tpms) <- rownames(ov_tpms_1)
for (i in 1:length(colnames(ov_tpms_1))) {ov_tpms[,i] <- NA}
colnames(ov_tpms) <- colnames(ov_tpms_1)
ov_tpms_1 <- as.data.frame(ov_tpms_1)
for (n in 1:length(colnames(ov_tpms_1))) {
  ov_tpms[,n] <- count_tpm_2(ov_tpms_1, exonic.gene.sizes)
  print(n)
}
thca_tpms_1 <- make_tcga_tpms("THCA")
thca_tpms <- as.data.frame(cbind(rownames(thca_tpms_1), rownames(thca_tpms_1)))
rownames(thca_tpms) <- rownames(thca_tpms_1)
for (i in 1:length(colnames(thca_tpms_1))) {thca_tpms[,i] <- NA}
colnames(thca_tpms) <- colnames(thca_tpms_1)
thca_tpms_1 <- as.data.frame(thca_tpms_1)
for (n in 1:length(colnames(thca_tpms_1))) {
  thca_tpms[,n] <- count_tpm_2(thca_tpms_1, exonic.gene.sizes)
  print(n)
}

scatter_table_oe_vs_ne_ig_tcga <- scatter_table_oe_vs_ne[c("ENSG00000124193.16", "ENSG00000116350", "ENSG00000251562",
                                                           "ENSG00000112081", "ENSG00000115875", "ENSG00000198625"),]
symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=strtrim(rownames(scatter_table_oe_vs_ne_ig_tcga), 15),mart= mart)
scatter_table_oe_vs_ne_ig_tcga <- scatter_table_oe_vs_ne_ig_tcga[!duplicated(rownames(scatter_table_oe_vs_ne_ig_tcga)),]
scatter_table_oe_vs_ne_ig_tcga$ensembl <- strtrim(rownames(scatter_table_oe_vs_ne_ig_tcga),15)
scatter_table_oe_vs_ne_ig_tcga$true <- scatter_table_oe_vs_ne_ig_tcga$ensembl%in%symbols$ensembl_gene_id
scatter_table_oe_vs_ne_ig_tcga <- scatter_table_oe_vs_ne_ig_tcga[!duplicated(scatter_table_oe_vs_ne_ig_tcga$ensembl),]
scatter_table_oe_vs_ne_ig_tcga <- scatter_table_oe_vs_ne_ig_tcga[order(rownames(scatter_table_oe_vs_ne_ig_tcga)),]
scatter_table_oe_vs_ne_ig_tcga$ensembl_2 <- symbols$ensembl_gene_id
scatter_table_oe_vs_ne_ig_tcga$symbol <- symbols$hgnc_symbol

ggplot(scatter_table_oe_vs_ne, aes(x=log2FoldChange_ne, y=log2FoldChange_oe)) +
  geom_pointdensity() +
  scale_color_viridis(option = "plasma", "neighbors") +
  theme_minimal() + geom_label_repel(scatter_table_oe_vs_ne_ig_tcga, 
                                     mapping = aes(label = symbol,
                                                   y = log2FoldChange_oe,
                                                   x = log2FoldChange_ne),
                                     box.padding = unit(0.35, "lines"),
                                     point.padding = unit(0.3, "lines"),
                                     force = 1, 
                                     fontface = "bold",
                                     size = 10/3)

dlbc_table <- make_table_for_scatters_tcga(dlbc_tpms)
gbm_table <- make_table_for_scatters_tcga(gbm_tpms)
lgg_table <- make_table_for_scatters_tcga(lgg_tpms)
thym_table <- make_table_for_scatters_tcga(thym_tpms)
kich_table <- make_table_for_scatters_tcga(kich_tpms)
ov_table <- make_table_for_scatters_tcga(ov_tpms)
thca_table <- make_table_for_scatters_tcga(thca_tpms)

kich_srsf4 <- plot_green_tcga(kich_table, kich_table$SRSF6, kich_table$SRSF4) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
kich_malat1 <- plot_green_tcga(kich_table, kich_table$SRSF6, kich_table$MALAT1) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
kich_srsf3 <- plot_green_tcga(kich_table, kich_table$SRSF6, kich_table$SRSF3) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
kich_srsf7 <- plot_green_tcga(kich_table, kich_table$SRSF6, kich_table$SRSF7) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
kich_chaf1a <- plot_green_tcga(kich_table, kich_table$SRSF6, kich_table$CHAF1A) + theme(axis.text.x=element_blank())
ov_srsf4 <- plot_green_tcga(ov_table, ov_table$SRSF6, ov_table$SRSF4) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
ov_malat1 <- plot_green_tcga(ov_table, ov_table$SRSF6, ov_table$MALAT1) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
ov_srsf3 <- plot_green_tcga(ov_table, ov_table$SRSF6, ov_table$SRSF3) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
ov_srsf7 <- plot_green_tcga(ov_table, ov_table$SRSF6, ov_table$SRSF7) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
ov_chaf1a <- plot_green_tcga(ov_table, ov_table$SRSF6, ov_table$CHAF1A) + theme(axis.text.x=element_blank())
thca_srsf4 <- plot_green_tcga(thca_table, thca_table$SRSF6, thca_table$SRSF4) + theme(axis.text.y=element_blank())
thca_malat1 <- plot_green_tcga(thca_table, thca_table$SRSF6, thca_table$MALAT1) + theme(axis.text.y=element_blank())
thca_srsf3 <- plot_green_tcga(thca_table, thca_table$SRSF6, thca_table$SRSF3) + theme(axis.text.y=element_blank())
thca_srsf7 <- plot_green_tcga(thca_table, thca_table$SRSF6, thca_table$SRSF7) + theme(axis.text.y=element_blank())
thca_chaf1a <- plot_green_tcga(thca_table, thca_table$SRSF6, thca_table$CHAF1A)

dlbc_srsf4 <- plot_red_tcga(dlbc_table, dlbc_table$SRSF6, dlbc_table$SRSF4) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
dlbc_malat1 <- plot_red_tcga(dlbc_table, dlbc_table$SRSF6, dlbc_table$MALAT1) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
dlbc_srsf3 <- plot_red_tcga(dlbc_table, dlbc_table$SRSF6, dlbc_table$SRSF3) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
dlbc_srsf7 <- plot_red_tcga(dlbc_table, dlbc_table$SRSF6, dlbc_table$SRSF7) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
dlbc_chaf1a <- plot_red_tcga(dlbc_table, dlbc_table$SRSF6, dlbc_table$CHAF1A) + theme(axis.text.x=element_blank())
gbm_srsf4 <- plot_red_tcga(gbm_table, gbm_table$SRSF6, gbm_table$SRSF4) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
gbm_malat1 <- plot_red_tcga(gbm_table, gbm_table$SRSF6, gbm_table$MALAT1) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
gbm_srsf3 <- plot_red_tcga(gbm_table, gbm_table$SRSF6, gbm_table$SRSF3) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
gbm_srsf7 <- plot_red_tcga(gbm_table, gbm_table$SRSF6, gbm_table$SRSF7) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
gbm_chaf1a <- plot_red_tcga(gbm_table, gbm_table$SRSF6, gbm_table$CHAF1A) + theme(axis.text.x=element_blank())
lgg_srsf4 <- plot_red_tcga(lgg_table, lgg_table$SRSF6, lgg_table$SRSF4) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
lgg_malat1 <- plot_red_tcga(lgg_table, lgg_table$SRSF6, lgg_table$MALAT1) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
lgg_srsf3 <- plot_red_tcga(lgg_table, lgg_table$SRSF6, lgg_table$SRSF3) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
lgg_srsf7 <- plot_red_tcga(lgg_table, lgg_table$SRSF6, lgg_table$SRSF7) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
lgg_chaf1a <- plot_red_tcga(lgg_table, lgg_table$SRSF6, lgg_table$CHAF1A) + theme(axis.text.x=element_blank())
thym_srsf4 <- plot_red_tcga(thym_table, thym_table$SRSF6, thym_table$SRSF4) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
thym_malat1 <- plot_red_tcga(thym_table, thym_table$SRSF6, thym_table$MALAT1) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
thym_srsf3 <- plot_red_tcga(thym_table, thym_table$SRSF6, thym_table$SRSF3) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
thym_srsf7 <- plot_red_tcga(thym_table, thym_table$SRSF6, thym_table$SRSF7) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
thym_chaf1a <- plot_red_tcga(thym_table, thym_table$SRSF6, thym_table$CHAF1A) + theme(axis.text.x=element_blank())

srsf4 <- textGrob("SRSF4", gp=gpar(fontface="bold"))
malat1 <- textGrob("MALAT1", gp=gpar(fontface="bold"))
srsf3 <- textGrob("SRSF3", gp=gpar(fontface="bold"))
srsf7 <- textGrob("SRSF7", gp=gpar(fontface="bold"))
chaf1a <- textGrob("MDM4", gp=gpar(fontface="bold"))
dlbc <- textGrob("DLBC", rot = 90, gp=gpar(fontface="bold"))
gbm <- textGrob("GBM", rot = 90, gp=gpar(fontface="bold"))
lgg <- textGrob("LGG", rot = 90, gp=gpar(fontface="bold"))
thym <- textGrob("THYM", rot = 90, gp=gpar(fontface="bold"))
kich <- textGrob("KICH", rot = 90, gp=gpar(fontface="bold"))
ov <- textGrob("OV", rot = 90, gp=gpar(fontface="bold"))
thca <- textGrob("THCA", rot = 90, gp=gpar(fontface="bold"))
blanc_text <- textGrob("", rot = 90, gp=gpar(fontface="bold"))

grid.arrange(blanc_text,srsf4, malat1,srsf3,srsf7,chaf1a,
             dlbc,dlbc_srsf4, dlbc_malat1,dlbc_srsf3,dlbc_srsf7,dlbc_chaf1a,
             gbm,gbm_srsf4,gbm_malat1,gbm_srsf3,gbm_srsf7,gbm_chaf1a,
             lgg,lgg_srsf4,lgg_malat1,lgg_srsf3,lgg_srsf7,lgg_chaf1a,
             thym,thym_srsf4,thym_malat1,thym_srsf3,thym_srsf7,thym_chaf1a,
             kich,kich_srsf4,kich_malat1,kich_srsf3,kich_srsf7,kich_chaf1a,
             ov,ov_srsf4,ov_malat1,ov_srsf3,ov_srsf7,ov_chaf1a,
             thca,thca_srsf4,thca_malat1,thca_srsf3,thca_srsf7,thca_chaf1a,
             nrow=8,
             widths = c(1,5,5,5,5,6),
             heights = c(1,5,5,5,5,5,5,6),
             as.table=TRUE)