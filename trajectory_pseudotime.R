library(Seurat)
library(reticulate)
library(slingshot)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(Matrix)
library(harmony)
library(SeuratWrappers)

setwd("/Users/hollywelfley/Desktop/New_Aspirate/")

#subset on myeloid of interest
aspirate=readRDS("Objects/cohort_labeltrans_singler.RDS")
types=unique(aspirate$FinalLabel)
keep=types[!types %in% c("Ciliated","Cycling","Basal","Neutrophil")]
myeloid=subset(aspirate, FinalLabel %in% keep)

#analysis workflow
myeloid <-FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 3000)
myeloid <- ScaleData(myeloid, verbose = FALSE)
myeloid <- RunPCA(myeloid, pc.genes = myeloid@var.genes, npcs = 20, verbose = FALSE)

myeloid <- myeloid %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(myeloid, 'harmony')

myeloid=RunUMAP(myeloid, reduction = "harmony", dims = 1:20)
myeloid=FindNeighbors(myeloid, reduction = "umap", dims = 1:2)
myeloid=FindClusters(myeloid, resolution = 0.09, algorithm = "4") 
myeloid=RunPHATE(myeloid, reduction="harmony")

saveRDS(myeloid, file="Objects/aspirate_myeloid_forpseudotime_030322.RDS")

#supplement plots
pdf(file="Supplement/aspiratemyeloid_umap.pdf", height=6, width=8)
DimPlot(myeloid, reduction = "umap", pt.size = 0.01)
dev.off()

pdf(file="Supplement/aspiratemyeloid_phate_bycluster.pdf", height=6, width=8)
DimPlot(myeloid, reduction = "phate", pt.size = 0.01)
dev.off()

# make object without cluster 4 (bridging cells)
notfour=subset(myeloid, seurat_clusters!=4)

#run slingshot for curve 1 and 2 
sds <- slingshot(Embeddings(notfour, "phate"), clusterLabels = notfour$merged, start.clus="3_6", end.clus=c(9,4))

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_clust <- cell_pal(notfour$merged, hue_pal())
cell_colors <- cell_pal(notfour$FinalLabel, c("#E0607E","#077187","#54BFB7","#3E5496","#8290BB"))

plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')

plot(reducedDim(sds), col = c("#E0607E","#077187","#54BFB7","#3E5496","#8290BB"), pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')

nc <- 2
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}
pdf(file="Figure_Drafts/aspirate_myeloid_pseudotimecurves_nocluster4.pdf",height=7,width=8)
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 2, col = 'black')
  scale_fill_viridis_c(option="magma")
}
dev.off()

curve=as.data.frame(pt)
colnames(curve)=c("Curve1","Curve2")

notfour=AddMetaData(notfour, curve)

counts=notfour@assays$RNA@counts
counts=counts[rowSums(counts)!=0,]
goodgenes=rownames(counts)
nozero=notfour[rowSums(notfour@assays$RNA@counts)!=0,]

saveRDS(nozero, "Objects/myeloid_nocluster4_withcurves_withoutzeros_030322.RDS")

c1=rownames(curve)[!is.na(curve$Curve1)]
c2=rownames(curve)[!is.na(curve$Curve2)]

curveone=nozero[,!colnames(nozero) %in% c1]
curvetwo=nozero[,!colnames(nozero) %in% c2]

curveone=subset(nozero, cells=c1)
curveone=curveone[rowSums(curveone@assays$RNA@counts)!=0,]

curvetwo=subset(nozero, cells=c2)
curvetwo=curvetwo[rowSums(curvetwo@assays$RNA@counts)!=0,]

saveRDS(curveone, file="Objects/myeloid_nocluster4_curveone_nozero_seuratobj.RDS")
saveRDS(curvetwo, file="Objects/myeloid_nocluster4_curvetwo_nozero_seuratobj.RDS")

saveRDS(notfour, file="Objects/myeloid_nocluster4_seuratobj.RDS")
saveRDS(sds, file="Objects/myeloid_nocluster4_slingshotobj.RDS")

########REMOVE MONO AND MAC#############################

mtypes=unique(myeloid$FinalLabel)
mkeep=mtypes[!mtypes %in% c("Monocyte","Macrophage")]
nomm=subset(myeloid, FinalLabel %in% mkeep)
DimPlot(nomm, reduction = "phate", group.by = "merged", label=T)

newsds <- slingshot(Embeddings(nomm, "phate"), clusterLabels = nomm$merged, start.clus="2_7", end.clus=4)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

newcell_colors <- cell_pal(nomm$FinalLabel, brewer_pal("qual", "Set2"))
newcell_colors_clust <- cell_pal(nomm$merged, hue_pal())

plot(reducedDim(newsds), col = newcell_colors, pch = 16, cex = 0.5)
lines(newsds, lwd = 2, type = 'lineages', col = 'black')

nc <- 2
newpt <- slingPseudotime(newsds)
newnms <- colnames(newpt)
newnr <- ceiling(length(newnms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(newnr, nc))
for (i in newnms) {
  colors <- pal[cut(newpt[,i], breaks = 100)]
  plot(reducedDim(newsds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(newsds, lwd = 2, col = 'black', type = 'lineages')
}
pdf(file="Figure_Drafts/aspirate_myeloid_pseudotimecurves_nomonocyte_nomacrophage.pdf",height=7,width=8)
for (i in newnms) {
  colors <- pal[cut(newpt[,i], breaks = 100)]
  plot(reducedDim(newsds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(newsds, lwd = 2, col = 'black')
  scale_fill_viridis_c(option="magma")
}
dev.off()

newcurve=as.data.frame(newpt)
colnames(newcurve)="Curve3"

nomm=AddMetaData(nomm, newcurve)

newcounts=nomm@assays$RNA@counts
newcounts=newcounts[rowSums(newcounts)!=0,]
newgoodgenes=rownames(newcounts)
newnozero=nomm[rowSums(nomm@assays$RNA@counts)!=0,]

saveRDS(newsds, file="Objects/myeloid_no_mono_mac_slingshotobj.RDS")

saveRDS(newnozero, "Objects/myeloid_nomonocyte_nomacrophage_withcurve3_withoutzeros_030322.RDS")

################################################
################################################
################################################
################################################
################################################
############## MONOCLE DEG ####################
################################################
################################################
################################################
################################################
################################################
################################################
library(Seurat)
library(monocle)
library(viridis)
library(Matrix)
one_expr_matrix=curveone@assays$RNA@counts
one_expr_matrix=as.matrix(one_expr_matrix)

one_meta_info=curveone@meta.data 
one_gene_annotations=as.data.frame(rownames(curveone@assays$RNA@counts))
colnames(one_gene_annotations)='gene_short_name'
rownames(one_gene_annotations)=rownames(curveone@assays$RNA@counts)
one_pd <- new("AnnotatedDataFrame", data = one_meta_info)
one_fd <- new("AnnotatedDataFrame", data = one_gene_annotations)
identical(colnames(one_expr_matrix), rownames(one_pd))
#TRUE
one_cds<- newCellDataSet(as.matrix(one_expr_matrix),
                         phenoData = one_pd, featureData = one_fd)

pData(one_cds)$Pseudotime=pData(one_cds)$Curve1
pData(one_cds)$SubjectID=pData(one_cds)$orig.ident
one_cds=estimateSizeFactors(one_cds)
one_cds=estimateDispersions(one_cds)

saveRDS(one_cds, "Objects/cds_curve1_030322.RDS")

#

two_expr_matrix=curvetwo@assays$RNA@counts
two_expr_matrix=as.matrix(two_expr_matrix)

two_meta_info=curvetwo@meta.data 
two_gene_annotations=as.data.frame(rownames(curvetwo@assays$RNA@counts))
colnames(two_gene_annotations)='gene_short_name'
rownames(two_gene_annotations)=rownames(curvetwo@assays$RNA@counts)
two_pd <- new("AnnotatedDataFrame", data = two_meta_info)
two_fd <- new("AnnotatedDataFrame", data = two_gene_annotations)
identical(colnames(two_expr_matrix), rownames(two_pd))
#TRUE
two_cds<- newCellDataSet(as.matrix(two_expr_matrix),
                         phenoData = two_pd, featureData = two_fd)

pData(two_cds)$Pseudotime=pData(two_cds)$Curve2
pData(two_cds)$SubjectID=pData(two_cds)$orig.ident
two_cds=estimateSizeFactors(two_cds)
two_cds=estimateDispersions(two_cds)

saveRDS(two_cds, "Objects/cds_curve2_030322.RDS")

#

three_expr_matrix=curvethree@assays$RNA@counts
three_expr_matrix=as.matrix(three_expr_matrix)

three_meta_info=curvethree@meta.data 
three_gene_annotations=as.data.frame(rownames(curvethree@assays$RNA@counts))
colnames(three_gene_annotations)='gene_short_name'
rownames(three_gene_annotations)=rownames(curvethree@assays$RNA@counts)
three_pd <- new("AnnotatedDataFrame", data = three_meta_info)
three_fd <- new("AnnotatedDataFrame", data = three_gene_annotations)
identical(colnames(three_expr_matrix), rownames(three_pd))
#TRUE
three_cds<- newCellDataSet(as.matrix(three_expr_matrix),
                           phenoData = three_pd, featureData = three_fd)

pData(three_cds)$Pseudotime=pData(three_cds)$Curve3
pData(three_cds)$SubjectID=pData(three_cds)$orig.ident
three_cds=estimateSizeFactors(three_cds)
three_cds=estimateDispersions(three_cds)

saveRDS(three_cds, "Objects/cds_curve3_030322.RDS")



#one with subject covariate
diff_test_covar1=differentialGeneTest(one_cds, fullModelFormulaStr = "~sm.ns(Pseudotime) + SubjectID",
                                      reducedModelFormulaStr = "~SubjectID", relative_expr = TRUE, cores = 1,
                                      verbose = FALSE)

#two with subject covariate
diff_test_covar2=differentialGeneTest(two_cds, fullModelFormulaStr = "~sm.ns(Pseudotime) + SubjectID",
                                      reducedModelFormulaStr = "~SubjectID", relative_expr = TRUE, cores = 1,
                                      verbose = FALSE)

#three with subject covariate
diff_test_covar3=differentialGeneTest(three_cds, fullModelFormulaStr = "~sm.ns(Pseudotime) + SubjectID",
                                      reducedModelFormulaStr = "~SubjectID", relative_expr = TRUE, cores = 1,
                                      verbose = FALSE)


diff_test_covar1=diff_test_covar1[!is.na(diff_test_covar1$qval),]
diff_test_covar1=diff_test_covar1[!is.na(diff_test_covar1$gene_short_name),]

diff_test_covar2=diff_test_covar2[!is.na(diff_test_covar2$qval),]
diff_test_covar2=diff_test_covar2[!is.na(diff_test_covar2$gene_short_name),]

diff_test_covar3=diff_test_covar3[!is.na(diff_test_covar3$qval),]
diff_test_covar3=diff_test_covar3[!is.na(diff_test_covar3$gene_short_name),]

write.table(diff_test_covar1,"trajectory_curve1_difftest_withcovar_nona_041222.txt",sep="\t",row.names=F,quote=F)
write.table(diff_test_covar2,"trajectory_curve2_difftest_withcovar_nona_041222.txt",sep="\t",row.names=F,quote=F)
write.table(diff_test_covar3,"trajectory_curve3_difftest_withcovar_nona_041222.txt",sep="\t",row.names=F,quote=F)

one_cosig_genes<- subset(diff_test_covar1, qval < 0.05)
two_cosig_genes<- subset(diff_test_covar2, qval < 0.05)
three_cosig_genes<- subset(diff_test_covar3, qval < 0.05)

one_sig_gene_names <- row.names(one_cosig_genes)
two_sig_gene_names <- row.names(two_cosig_genes)
three_sig_gene_names <- row.names(three_cosig_genes)

one_cds_subset=one_cds[one_cds@featureData@data$gene_short_name %in% one_sig_gene_names,]
saveRDS(one_cds_subset,"Objects/cds_subset_sigDEG_curve1_030322.RDS")

two_cds_subset=two_cds[two_cds@featureData@data$gene_short_name %in% two_sig_gene_names,]
saveRDS(two_cds_subset,"Objects/cds_subset_sigDEG_curve2_030322.RDS")

three_cds_subset=three_cds[three_cds@featureData@data$gene_short_name %in% three_sig_gene_names,]
saveRDS(three_cds_subset,"Objects/cds_subset_sigDEG_curve3_030322.RDS")
three_cds_subset=readRDS("Objects/cds_subset_sigDEG_curve3_030322.RDS")

#
##################################
#heatmap
##################################

pdf(file="pseudotime_heatmap_curve1_030322.pdf", height=6, width=6)
plot_pseudotime_heatmap(one_cds_subset,num_clusters = 4,cores = 1,show_rownames = F, hmcols=viridis(200,option="magma"))
dev.off()

pdf(file="pseudotime_heatmap_curve2_030322.pdf", height=6, width=6)
plot_pseudotime_heatmap(two_cds_subset,num_clusters = 3,cores = 1,show_rownames = F, hmcols=viridis(200,option="magma"))
dev.off()

pdf(file="pseudotime_heatmap_curve3_030322.pdf", height=6, width=6)
plot_pseudotime_heatmap(three_cds_subset,num_clusters = 3,cores = 1,show_rownames = F, hmcols=viridis(200,option="magma"))
dev.off()

#
one_heatmap=plot_pseudotime_heatmap(one_cds_subset,num_clusters = 4,cores = 1,show_rownames = F, return_heatmap=T)
two_heatmap=plot_pseudotime_heatmap(two_cds_subset,num_clusters = 3,cores = 1,show_rownames = F, return_heatmap=T)
three_heatmap=plot_pseudotime_heatmap(three_cds_subset,num_clusters = 3,cores = 1,show_rownames = F, return_heatmap=T)

one_annotation_row <- data.frame(Cluster=factor(cutree(one_heatmap$tree_row, 4)))
one_annotation_row$genename=rownames(one_annotation_row)
one_data=merge(one_cosig_genes, one_annotation_row, by.x="gene_short_name", by.y="genename", all.x=T, sort=F)
write.table(one_data, file="pseudotime_curve1_diffexpgenes_4clust_withclusterinfo_030322.txt", row.names=T, col.names = T, quote=T, sep="\t")

two_annotation_row <- data.frame(Cluster=factor(cutree(two_heatmap$tree_row, 3)))
two_annotation_row$genename=rownames(two_annotation_row)
two_data=merge(two_cosig_genes, two_annotation_row, by.x="gene_short_name", by.y="genename", all.x=T, sort=F)
write.table(two_data, file="pseudotime_curve2_diffexpgenes_3clust_withclusterinfo_030322.txt", row.names=T, col.names = T, quote=T, sep="\t")

three_annotation_row <- data.frame(Cluster=factor(cutree(three_heatmap$tree_row, 3)))
three_annotation_row$genename=rownames(three_annotation_row)
three_data=merge(three_cosig_genes, three_annotation_row, by.x="gene_short_name", by.y="genename", all.x=T, sort=F)
write.table(three_data, file="pseudotime_curve3_diffexpgenes_3clust_withclusterinfo_030322.txt", row.names=T, col.names = T, quote=T, sep="\t")
three_data=read.table(file="pseudotime_curve3_diffexpgenes_3clust_withclusterinfo_030322.txt", header=T)

onegene=data.frame(Gene=one_data$gene_short_name,
                   Pval=one_data$pval,
                   Qval=one_data$qval,
                   Cluster=one_data$Label)
write.table(onegene,"Supplement/Tables/TableS7_am_trajectory_genes.txt",sep="\t", col.names = T,row.names=F,quote=F)

twogene=data.frame(Gene=two_data$gene_short_name,
                   Pval=two_data$pval,
                   Qval=two_data$qval,
                   Cluster=two_data$Label)
write.table(twogene,"Supplement/Tables/TableS9_im_trajectory_genes.txt",sep="\t", col.names = T,row.names=F,quote=F)

threegene=data.frame(Gene=three_data$gene_short_name,
                     Pval=three_data$pval,
                     Qval=three_data$qval,
                     Cluster=three_data$Label)
write.table(threegene,"Supplement/Tables/TableSX_bridging_trajectory_genes.txt",sep="\t", col.names = T,row.names=F,quote=F)

########
#get gene lists for strict (qval<0.05) and relaxed (pval<0.05)
########
amstrict=diff_test_covar1$gene_short_name[diff_test_covar1$qval<0.05]
#4240
amrelaxed=diff_test_covar1$gene_short_name[diff_test_covar1$pval<0.05]
#5324
imstrict=diff_test_covar2$gene_short_name[diff_test_covar2$qval<0.05]
#3698
imrelaxed=diff_test_covar2$gene_short_name[diff_test_covar2$pval<0.05]
#4921
bridgestrict=diff_test_covar3$gene_short_name[diff_test_covar3$qval<0.05]
#4567
bridgerelaxed=diff_test_covar3$gene_short_name[diff_test_covar3$pval<0.05]
#5599

### Generate gene lists of pval<0.05 for trajectory combos
amimrelaxed=union(amrelaxed, imrelaxed)
#7626
ambridgerelaxed=union(amrelaxed, bridgerelaxed)
#7041
bridgeimrelaxed=union(bridgerelaxed, imrelaxed)
#7610

#### Find trajectory specific genes (needs to be qval<0.05 and NOT pval<0.05 in the other two categories)
allbridgespecific=setdiff(x=bridgestrict, y=amimrelaxed)
#285
allamspecific=setdiff(x=amstrict, y=bridgeimrelaxed)
#259
allimspecific=setdiff(x=imstrict, y=ambridgerelaxed)
#666

############### filter allrelaxed genes by needing to have qval<0.05 in at least one category

amshared=intersect(amstrict, bridgeimrelaxed)
#3981
imshared=intersect(imstrict, ambridgerelaxed)
#3032
bridgeshared=intersect(bridgestrict, amimrelaxed)
#4282

##### total shared 

two_1=union(amshared, imshared)
#5188
two=union(two_1, bridgeshared)
#5451

#totals
amim=union(amstrict, imstrict)
#6113
three=union(bridgestrict, amim)
#6661

#######################
#make data frames
allamspecificdf=data.frame(gene_short_name=allamspecific, Label="AM_Specific")
allimspecificdf=data.frame(gene_short_name=allimspecific, Label="IM_Specific")
allbridgespecificdf=data.frame(gene_short_name=allbridgespecific, Label="Bridge_Specific")
amshareddf=data.frame(gene_short_name=amshared, Label="AM_Shared")
imshareddf=data.frame(gene_short_name=imshared, Label="IM_Shared")
bridgeshareddf=data.frame(gene_short_name=bridgeshared, Label="Bridge_Shared")
allshareddf=data.frame(gene_short_name=sharedanytwo, Label="Shared")

alldata=rbind(allamspecificdf, allimspecificdf, allbridgespecificdf, allshareddf)
formatalldata=alldata
colnames(formatalldata)=c("Gene","Label")
write.table(formatalldata,"pseudotime_enrichr_moredbs/RELAXED_AMIMBridgeSpecific_SharedAllThree_genes.txt",sep="\t",row.names=F,quote=F)

###################
###################
###################
###################
# GO Terms
###################
###################
###################
###################

dbs <- c("GO_Biological_Process_2021","WikiPathway_2021_Human",
         "Reactome_2016")

enricher = function(de_list){
  clusters = unique(de_list$Label)
  enrichout=data.frame(Term=character(),
                       Overlap=character(),
                       P.value=character(),
                       Adjusted.P.value=character(),
                       Old.P.value=character(),
                       Old.Adjusted.P.value=character(),
                       Odds.Ratio=character(),
                       Combined.Score=character(),
                       Genes=character(),
                       stringsAsFactors=FALSE)
  
  for(cluster in clusters){
    gex.subset = de_list[de_list$Label == cluster,]
    if(nrow(gex.subset) < 1){next}
    curr_enriched = enrichr(gex.subset$gene_short_name, dbs)
    for(j in 1:length(curr_enriched)){
      curr_db = curr_enriched[[j]]
      if(nrow(curr_db) < 1){next}
      curr_db$Cluster = cluster
      curr_db$database = names(curr_enriched)[j]
      enrichout = rbind(enrichout,curr_db)
    }
  }
  enrichout$gene_hits = as.numeric(gsub("/.*","", enrichout$Overlap))
  enrichout$pathway_genes = as.numeric(gsub(".*/","", enrichout$Overlap))
  enrichout_filtered=enrichout[enrichout$Adjusted.P.value<0.05,]
  enrichout_filtered=enrichout_filtered[enrichout_filtered$gene_hits>1,]
  enrichout_formatted = enrichout_filtered[order(enrichout_filtered$Adjusted.P.value),c(10,11,1,12,13,3,4,7,8,9)]
  return(enrichout_formatted)
}

alldata_paths = enricher(alldata)
write.table(alldata_paths,"pseudotime_enrichr_moredbs/RELAXED_AMvsIMvsBridge_all_pathways.txt",sep="\t",row.names=F,quote=F)

AM_paths = enricher(one_data)
write.table(AM_paths,"pseudotime_enrichr_moredbs/AM_pathways.txt",sep="\t",row.names=F,quote=F)
IM_paths = enricher(two_data)
write.table(IM_paths,"pseudotime_enrichr_moredbs/IM_pathways.txt",sep="\t",row.names=F,quote=F)
bridge_paths = enricher(three_data)
write.table(bridge_paths,"pseudotime_enrichr_moredbs/Bridging_pathways.txt",sep="\t",row.names=F,quote=F)
cluster4genes_moreinfo=read.table("bridgingspecificgenes_cluster4_bin19to62_pqval_location.txt",header=T)
forbridge=data.frame(gene_short_name=cluster4genes_moreinfo$Gene, Label="BridgeSpecific")
specificbridge_paths = enricher(forbridge)
write.table(specificbridge_paths,"pseudotime_enrichr_moredbs/Bridging_cluster4_bin19to62_pathways.txt",sep="\t",row.names=F,quote=F)

#####
############
#################get bridging specific and cluster 4 info
############
#####

bridgelocation=read.table(file="Supplement/Tables/TableSX_bridging_trajectory_genes.txt", header=T)
bridgespecificwithlocation=bridgelocation[bridgelocation$Gene %in% bridgespecific,]
write.table(bridgespecificwithlocation,"pseudotime_enrichr_moredbs/RELAXED_BridgeSpecific_genes_withlocation.txt",sep="\t",row.names=F,quote=F)

bridgespecificdata=read.table(file="pseudotime_enrichr_moredbs/RELAXED_BridgeSpecific_genes_withlocation.txt", header=T)
bridgespecificgenes=bridgespecificdata$Gene
bridgesubset=three_cds_subset[three_cds_subset@featureData@data$gene_short_name %in% bridgespecificgenes,]
pdf(file="pseudotime_heatmap_bridgespecific.pdf", height=6, width=6)
plot_pseudotime_heatmap(bridgesubset,num_clusters = 3,cores = 1,show_rownames = F, hmcols=viridis(200,option="magma"))
dev.off()

bridge_heatmap=plot_pseudotime_heatmap(bridgesubset,num_clusters = 3,cores = 1,show_rownames = F, return_heatmap=T)

newdata <- data.frame(Pseudotime = seq(min(pData(bridgesubset)$Pseudotime), max(pData(bridgesubset)$Pseudotime),length.out = 100)) 
m <- genSmoothCurves(bridgesubset, cores=1, trend_formula = '~sm.ns(Pseudotime, df=3)',  
                     relative_expr = T, new_data = newdata)

maxer = function(x){
  maxnow = max(x)
  maxes = which(x == maxnow)
  if(length(maxes) == 1){
    return(maxes)
  }else{
    maxes = paste(maxes,collapse="_")
  }
}

maxed = apply(m,1,maxer)
#cluster 4 pseudotime range = bins 19:62
rownumsgenes = names(which(maxed > 19 & maxed < 62))

hm=as.data.frame(t(m))
hm$Cluster=rownames(hm)

hmaxed=as.data.frame(maxed)
hmaxed$Gene=rownames(hmaxed)
rownames(hmaxed)=NULL
colnames(hmaxed)=c("Cluster","Gene")

withexp=data.frame("Gene"=as.character(), "Cluster"=as.character(), "ExpValue"==as.character(), stringsAsFactors = F)

for(name in hmaxed$Gene){
  maxbin=hmaxed$Cluster[hmaxed$Gene==name]
  thishm=hm[,which(colnames(hm)==name)]
  keep=thishm[maxbin]
  keepme=data.frame("Gene"=name, "Cluster"=maxbin, "ExpValue"=keep)
  withexp=rbind(withexp, keepme)
}

identical(hmaxed, withexp[,c(2,1)])
#TRUE

colnames(withexp)=c("Gene","MaxSpline_Bin","MaxSpline_Expression")
withexp$BridgingRegion=withexp$MaxSpline_Bin
good=which(withexp$BridgingRegion> 19 & withexp$BridgingRegion < 62)
withexp$BridgingRegion[good]="Yes"
withexp$BridgingRegion[withexp$BridgingRegion!="Yes"]="No"

write.table(withexp, file="/Users/hollywelfley/Desktop/New_Aspirate/pseudotime_enrichr_moredbs/bridgingspecificgenes_withbinandexpinfo_forssix.txt",row.names = F, col.names = T, quote=T, sep="\t")
write.table(rownumsgenes, file="bridgingspecificgenes_cluster4_bin19to62.txt",row.names = F, col.names = F, quote=T, sep="\t")

