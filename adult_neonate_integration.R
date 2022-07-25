library(Seurat)
library(harmony)
library(leiden)
library(Libra)
library(gplots)
library(RColorBrewer)
setwd("/xdisk/darrenc/hwelfley/redo/integration")

#read in datasets
tgen=readRDS("/xdisk/darrenc/hwelfley/aspirate_project/lung_ref/references/tgen_controllung.RDS")
hlca=readRDS("/xdisk/darrenc/hwelfley/aspirate_project/lung_ref/references/lung_counts_ref_sparse_allPatients_RunUMAP_061020.rds")
aspirate=readRDS("/xdisk/darrenc/hwelfley/redo/labeltransfer/cohort_labeltrans_singler.RDS")

#subset on genes in all datasets
av=as.vector(aspirate@assays$RNA@counts@Dimnames[[1]])
hv=as.vector(hlca@assays$RNA@counts@Dimnames[[1]])
tv=as.vector(tgen@assays$RNA@counts@Dimnames[[1]])

lunggene=intersect(tv, hv)
igenes=intersect(lunggene, av)
length(igenes)

it=subset(tgen, features=igenes)
ih=subset(hlca, features=igenes)
ia=subset(aspirate, features=igenes)

#carry over cell type and dataset info to shared metadata columns
#make large categories for later subsetting
it$CombinedLabels=it$celltype
it$SampleType="TGEN_Lung"
it$GeneralType=it$population
it$FilterMyeloid=it$CombinedLabels
it$FilterMyeloid[it$GeneralType=="Epithelial"]="Not_Myeloid"
it$FilterMyeloid[it$GeneralType=="Mesenchymal"]="Not_Myeloid"
it$FilterMyeloid[it$GeneralType=="Endothelial"]="Not_Myeloid"
it$FilterMyeloid[it$FilterMyeloid=="Proliferating Macrophages"]="Proliferating"
it$FilterMyeloid[it$FilterMyeloid %in% c("Proliferating T Cells","NK Cells","Mast Cells","Plasma Cells","B Cells","T Cells")]="Not_Myeloid"
it$FilterMyeloid[it$FilterMyeloid %in% c("pDCs","cDCs")]="Dendritic" 
it$FilterMyeloid[it$FilterMyeloid %in% c("Monocytes","Macrophages")]="Myeloid" 

ih$CombinedLabels=ih$free_annotation
ih$SampleType="HLCA_Lung"
ih$GeneralType=as.vector(ih$compartment)
ih$FilterMyeloid=as.vector(ih$CombinedLabels)
ih$FilterMyeloid[ih$GeneralType=="Epithelial"]="Not_Myeloid"
ih$FilterMyeloid[ih$GeneralType=="Stromal"]="Not_Myeloid"
ih$FilterMyeloid[ih$GeneralType=="Endothelial"]="Not_Myeloid"
ih$FilterMyeloid[ih$FilterMyeloid %in% c("IGSF21+ Dendritic","Myeloid Dendritic Type 1","Myeloid Dendritic Type 2", "Plasmacytoid Dendritic","EREG+ Dendritic",
                                         "TREM2+ Dendritic")]="Dendritic"
ih$FilterMyeloid[ih$FilterMyeloid %in% c("B","CD8+ Naive T","CD4+ Naive T","CD4+ Memory/Effector T","CD8+ Memory/Effector T","Proliferating NK/T","Natural Killer T","Natural Killer","Basophil/Mast 1","Basophil/Mast 2", "Platelet/Megakaryocyte","Plasma")]="Not_Myeloid"
ih$FilterMyeloid[ih$FilterMyeloid=="Proliferating Macrophage"]="Proliferating"
ih$FilterMyeloid[!ih$FilterMyeloid %in% c("Proliferating","Dendritic","Not_Myeloid")]="Myeloid"

ia$CombinedLabels=ia$FinalLabel
ia$SampleType="Aspirate"
ia$FilterMyeloid=ia$CombinedLabels
ia$FilterMyeloid[ia$FilterMyeloid %in% c("Basal","Ciliated","Neutrophil")]="Not_Myeloid"
ia$FilterMyeloid[ia$FilterMyeloid=="Cycling"]="Proliferating"
ia$FilterMyeloid[!ia$FilterMyeloid %in% c("Proliferating","Not_Myeloid")]="Myeloid"
ia$GeneralType=ia$FilterMyeloid
ia$GeneralType[ia$FinalLabel %in% c("Basal","Ciliated")]="Epithelial"
ia$GeneralType[ia$FinalLabel=="Neutrophil"]="Granulocyte"

#merge and save
all=merge(it, y=c(ih, ia))
saveRDS(all, "new/aspirate_tgen_hlca_onlymerged.RDS")  

# all cell types umap
all <- NormalizeData(all,verbose = FALSE)
all <-FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(all, verbose = FALSE)
all <- RunPCA(all, pc.genes = all@var.genes, npcs = 30, verbose = FALSE)

all <- all %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony=100)
harmony_embeddings <- Embeddings(all, 'harmony')

all <- all %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4, algorithm = "4") %>% 
  identity()

pdf(file="new/aspirate_lungs_umaps.pdf", height=12, width=40)
h1=DimPlot(all, reduction="umap", label=T, pt.size = 0.01)
h2=DimPlot(all, reduction="umap", label=T, group.by = "SampleType", pt.size = 0.01)
h3=DimPlot(all, reduction="umap", label=T, group.by = "CombinedLabels", pt.size = 0.01)+NoLegend()
h1+h2+h3
dev.off()

saveRDS(all, "new/aspirate_tgen_hlca_harmony.RDS")  

####################################
####################################
### Subset on myeloid cells of interest
####################################
####################################

all$BetterLabel=all$CombinedLabels
all$BetterLabel[all$BetterLabel %in% c("Monocytes","OLR1+ Classical Monocyte","Classical Monocyte")]="Monocyte"
all$BetterLabel[all$BetterLabel=="Macrophages"]="Macrophage"
goodcells=as.vector(unique(all$BetterLabel[all$FilterMyeloid=="Myeloid"]))
#"Monocyte"                   "Macrophage"                
#"Intermediate Monocyte"      "Nonclassical Monocyte"     
#"Intermediate I (VDAC1)"     "Intermediate II (MIR222HG)"
#"Intermediate III (BEST1)"
myeloid=subset(all, BetterLabel %in% goodcells)

#subset to same number
table(myeloid$SampleType)
#Aspirate HLCA_Lung TGEN_Lung 
#5161     18010     14858 

aspcell=rownames(myeloid@meta.data)[myeloid$SampleType=="Aspirate"]
set.seed(456)
hlcacell=sample(rownames(myeloid@meta.data)[myeloid$SampleType=="HLCA_Lung"], 5161)
tgencell=sample(rownames(myeloid@meta.data)[myeloid$SampleType=="TGEN_Lung"], 5161)

mono=subset(myeloid, cells=c(aspcell, hlcacell, tgencell))

# run analysis
mono <-FindVariableFeatures(mono, selection.method = "vst", nfeatures = 3000)
mono <- ScaleData(mono, verbose = FALSE)
mono <- RunPCA(mono, pc.genes = mono@var.genes, npcs = 20, verbose = FALSE)

mono <- mono %>% 
  RunHarmony(c("orig.ident","SampleType"), plot_convergence = TRUE, max.iter.harmony=100)
harmony_embeddings <- Embeddings(mono, 'harmony')

mono <- mono %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution=0.45, algorithm = "4") %>% 
  identity()

#relabel the clusters with new labels from analysis
mono$FigureLabel=as.vector(mono$seurat_clusters)
mono$FigureLabel[mono$FigureLabel %in% c(1,5,7,8,10,12)]="Macrophage"
mono$FigureLabel[mono$FigureLabel==4]="InterI_Mac"
mono$FigureLabel[mono$FigureLabel==9]="InterII_Mac"
mono$FigureLabel[mono$FigureLabel==2]="Neonate Inter"
mono$FigureLabel[mono$FigureLabel %in% c(3,11,6)]="Monocyte"

saveRDS(mono, file = "Objects/asp_hlca_tgen_subsam_031322.RDS")

#generate confusion matrix
dist=data.frame(unclass(table(mono$BetterLabel, mono$seurat_clusters)))
labdist=as.matrix(dist)
labsums=colSums(labdist)  
newlabdist=t(labdist)
labper=newlabdist/labsums
rownames(labper)=c(1:12)
newlabper=t(labper)
reordered=newlabper[,c(1,5,7,8,10,12, 3,6,11,2,4,9)]

pdf(file="supplement/adultvneolabeldist_acrosstotalcluster.pdf",height=8, width=10)
pheatmap(reordered,cluster_rows = F, cluster_cols = F,
         color=viridis(100,option="viridis"),
         border_color =NA,
         angle_col = 45,cellwidth = 20, cellheight = 20,
         show_rownames = T, show_colnames = T)
dev.off()

labpercent=labper*100
write.csv(labpercent, file="labeldist_acrosstotalcluster.csv", col.names = T, row.names = T)

#write figures
lesscolor=c("#BF40BF", "#00D1D1", "#3E5496", "#8290BB", "#8E2043")

pdf(file="Figure_Drafts/adult_neonate_umap_031322.pdf", height=6, width=8)
DimPlot(mono, reduction = "umap", group.by = "FigureLabel", pt.size = 0.01, cols = lesscolor)
dev.off()

pdf(file="Figure_Drafts/adult_neonate_celltypedist_barplot_031322.pdf", height=4, width=4)
barplot(plotper, col=lesscolor)
dev.off()

#######################################
#### Libra DEGs
#######################################

mono$replicate=mono$orig.ident
mono$label=mono$FigureLabel
mono$cell_type=mono$FigureLabel

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

clustermarkerer(mono,"vs_adultneomyeloid")

little=subset(mono, FigureLabel %in% c("InterI_Mac","InterII_Mac"))
clustermarkerer(little,"vs_intermac_eachother")

newonemac=read.table("Libra_Outs/InterI_Mac_vs_intermac_eachother_DE.txt", header=T)
supmacs=newonemac
supmacs$cell_type[supmacs$avg_logFC>0]="Mac_InterI"
supmacs$cell_type[supmacs$avg_logFC<0]="Mac_InterII"
supmacs$avg_logFC=abs(supmacs$avg_logFC)
supmacs=supmacs[supmacs$p_val_adj<0.05,]
bettersupmacs=data.frame(Gene=supmacs$gene,
                         AvgLogFC=supmacs$avg_logFC,
                         Pval=supmacs$p_val,
                         AdjPval=supmacs$p_val_adj,
                         CellType=supmacs$cell_type)
write.table(bettersupmacs,"Supplement/Tables/TableSx_maciI_vs_maciII_markers.txt",sep="\t", col.names = T,row.names=F,quote=F)
