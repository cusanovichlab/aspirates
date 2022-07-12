library(Seurat)
library(harmony)
library(gplots)
library(RColorBrewer)

setwd("/Users/hollywelfley/Desktop/New_Aspirate/")
cohort=readRDS("/Users/hollywelfley/Desktop/New_Aspirate/cohort_onlymerged_seuratobj_soupx_scrublet200minfeatmt20.RDS")

cohort <- Seurat::NormalizeData(cohort,verbose = FALSE)
cohort <-FindVariableFeatures(cohort, selection.method = "vst", nfeatures = 3000)
cohort <- ScaleData(cohort, verbose = FALSE)
cohort <- RunPCA(cohort, pc.genes = cohort@var.genes, npcs = 20, verbose = FALSE)

cohort <- RunHarmony(cohort, "orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(cohort, 'harmony')

cohort <- cohort %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "umap", dims = 1:2) %>% 
  FindClusters(resolution = 0.12, algorithm = "4") %>% 
  identity()

DimPlot(cohort, reduction="umap", label=T, pt.size=0.001)

saveRDS(cohort, file="cohort_3kvarfeat_20pc_2dim_res0p12.RDS")

###########
#Spearman correlation of clusters
###########

aves = AverageExpression(cohort,assays="RNA",group.by = "RNA_snn_res.0.12")
test_scale = scale(log1p(aves[[1]]))
aves_na = aves[[1]]
aves_na[aves_na == 0] = NA

aves.myeloid = aves[[1]][,c(1:7,9)]
aves.myeloid_na = aves.myeloid
aves.myeloid_na[aves.myeloid_na == 0] = NA
pca.cluster = prcomp(t(aves.myeloid[rowSums(aves.myeloid)>0,]),scale=T)


pdf(file="Recent_Figures/cohort_3k20pc_2dim0p12res_heatmap_clusters.pdf", height=8, width=8)
heatmap.2(cor(aves_na,method="spearman",use="pairwise.complete.obs"), col = brewer.pal(9,"Blues"), margins = c(14, 14), 
          density.info = "none", lhei=c(2, 8), trace= "none", #ColSideColors = assaycols,
          distfun=function(x) as.dist(1-abs(x)),
          hclustfun=function(x) hclust(x,method="ward.D2"),main="Ward clusters")
dev.off()

pdf(file="Recent_Figures/cohort_3k20pc_2dim0p12res_heatmap_myeloidonly_clusters.pdf", height=8, width=8)
heatmap.2(cor(aves.myeloid_na,method="spearman",use="pairwise.complete.obs"), col = brewer.pal(9,"Blues"), margins = c(14, 14), 
          density.info = "none", lhei=c(2, 8), trace= "none", #ColSideColors = assaycols,
          distfun=function(x) as.dist(1-abs(x)),
          hclustfun=function(x) hclust(x,method="ward.D2"),main="Ward clusters")
dev.off()

######
#merge clusters by correlation
######

cohort$merged=as.vector(cohort$seurat_clusters)
cohort$merged[cohort$merged==3]="3_6"
cohort$merged[cohort$merged==6]="3_6"
cohort$merged[cohort$merged==2]="2_7"
cohort$merged[cohort$merged==7]="2_7"

pdf(file="aspirate_umap_harmony_orig.ident.pdf", height=6, width=8)
DimPlot(cohort, reduction="umap", group.by = "orig.ident", pt.size=0.001)
dev.off()

saveRDS(cohort, file="cohort_3kvarfeat_20pc_2dim_res0p12_mergedclusters.RDS")

###########
#Libra differential expression
###########

library(Libra)

cohort$replicate=cohort$orig.ident
cohort$label=cohort$merged
cohort$cell_type=cohort$merged

clustermarkerer = function(seurat_obj,namer){
  meta_save = seurat_obj@meta.data
  cluster_combos = unique(meta_save$label)
  cell_types = meta_save[match(cluster_combos,meta_save$label),c("cell_type","label")]
  seurat_obj@meta.data$cell_type = "Markers"
  table_report = matrix(,,5)
  for(i in 1:length(cluster_combos)){
    print(cell_types[i,1])
    seurat_obj@meta.data$label = (meta_save$label != cluster_combos[i])+0
    temp_subset = subset(seurat_obj,cells=rownames(seurat_obj@meta.data)[meta_save$label == cluster_combos[i]])
    temp_gene_counts = rowSums(temp_subset@assays$RNA@counts>0)
    gene.list = names(temp_gene_counts[which(temp_gene_counts >= round(length(which(meta_save$label == cluster_combos[i]))*.05))])
    temp_subset2 = subset(seurat_obj,cells=rownames(seurat_obj@meta.data)[meta_save$label != cluster_combos[i]])
    temp_gene_counts2 = rowSums(temp_subset2@assays$RNA@counts>0)
    gene.list2 = names(temp_gene_counts2[which(temp_gene_counts2 >= round(length(which(meta_save$label != cluster_combos[i]))*.05))])
    final.gene.list = union(gene.list,gene.list2)
    seurat_obj = subset(seurat_obj,features = final.gene.list)
    DE = run_de(seurat_obj)
    DE_o = DE[order(DE$p_val_adj),]
    DE_o$cell_type = cell_types[i,1]
    DE_sig = DE_o[DE_o$p_val_adj < 0.05,]
    DE_pos = DE_sig[DE_sig$avg_logFC > 0,]
    print(paste0(nrow(DE_pos)," sig. markers"))
    print(paste0(nrow(DE_sig)-nrow(DE_pos)," sig. downregulated"))
    print(paste0("New marker is: ",DE_pos$gene[1]))
    write.table(DE_o,paste0("Libra_Outs/",cell_types[i,1],"_",namer,"_DE.txt"),row.names=F,quote=F,sep="\t")
    table_report = rbind(table_report,c(cell_types[i,1],nrow(DE_pos),nrow(DE_sig)-nrow(DE_pos),nrow(DE_o),DE_pos$gene[1]))
  }
  colnames(table_report) = c("Cluster","DE_up","DE_down","GenesTested","TopMarker")
  table_report = table_report[-1,]
  write.table(table_report,paste0("Libra_Outs/one_",namer,"_DE_table_report.txt"),row.names=F,quote=F,sep="\t")
}

clustermarkerer(cohort,"vs_All")

#see other methods of cell annotation (label transfer, SingleR) that were used to define cell labels
cohort$FinalLabel=cohort$merged
#ordered by size
cohort$FinalLabel[cohort$FinalLabel=="3_6"]="Monocyte"
cohort$FinalLabel[cohort$FinalLabel=="2_7"]="Intermediate I (VDAC1)"
cohort$FinalLabel[cohort$FinalLabel==1]="Neutrophil"
cohort$FinalLabel[cohort$FinalLabel==4]="Intermediate II (MIR222HG)"
cohort$FinalLabel[cohort$FinalLabel==5]="Intermediate III (BEST1)"
cohort$FinalLabel[cohort$FinalLabel==8]="Basal"
cohort$FinalLabel[cohort$FinalLabel==9]="Macrophage"
cohort$FinalLabel[cohort$FinalLabel==10]="Ciliated"
cohort$FinalLabel[cohort$FinalLabel==11]="Cycling"

DimPlot(cohort, reduction="umap", pt.size=0.001, group.by = "FinalLabel")

saveRDS(cohort, "cohort_3kvarfeat_20pc_2dim_res0p12_finallabels.RDS")

#merge DEGs/markers for supplement file

final=data.frame(cell_type=character(),
                     gene=character(),
                     avg_logFC=character(),
                     p_val=character(),
                     p_val_adj=character(),
                     de_family=character(),
                     de_method=character(),
                     de_type=character(),
                     stringsAsFactors=FALSE)

setwd("/Users/hollywelfley/Desktop/New_Aspirate/Libra_Outs/")
for (f in c("1_vs_All_DE.txt", "2_7_vs_All_DE.txt", "3_6_vs_All_DE.txt",
               "4_vs_All_DE.txt","5_vs_All_DE.txt","8_vs_All_DE.txt",
               "9_vs_All_DE.txt","10_vs_All_DE.txt","11_vs_All_DE.txt")){
  current=read.table(file=paste0(f), header = T)
  final=rbind(final, current)}
  
final=final[final$p_val_adj<0.05,]
final$cell_type[final$cell_type=="3_6"]="Monocyte"
final$cell_type[final$cell_type=="2_7"]="Intermediate I (VDAC1)"
final$cell_type[final$cell_type==1]="Neutrophil"
final$cell_type[final$cell_type==4]="Intermediate II (MIR222HG)"
final$cell_type[final$cell_type==5]="Intermediate III (BEST1)"
final$cell_type[final$cell_type==8]="Basal"
final$cell_type[final$cell_type==9]="Macrophage"
final$cell_type[final$cell_type==10]="Ciliated"
final$cell_type[final$cell_type==11]="Cycling"
  
supment=data.frame(Gene=final$gene,
                   AvgLogFC=final$avg_logFC,
                   Pval=final$p_val,
                   AdjPval=final$p_val_adj,
                   CellType=final$cell_type)

write.table(supment,"Supplement/Tables/TableS3_aspiratecells_degs.txt",sep="\t", col.names = T,row.names=F,quote=F)

  
