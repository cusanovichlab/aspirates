library(Seurat)
library(pheatmap)
library(viridis)

setwd("/xdisk/darrenc/hwelfley/redo/labeltransfer")
cohort=readRDS("cohort_3kvarfeat_20pc_2dim_res0p12_finallabels.RDS")

#Travaglini HLCA

hlca.ref.rna <-readRDS("/xdisk/darrenc/hwelfley/aspirate_project/lung_ref/references/lung_counts_ref_sparse_allPatients_RunUMAP_061020.rds")

# hlca label transfer
hlca.anchors <- FindTransferAnchors(reference = hlca.ref.rna, query = cohort, dims = 1:30)
hlca_free_annotation_predictions <- TransferData(anchorset = hlca.anchors, refdata = hlca.ref.rna$free_annotation, dims = 1:30)
colnames(hlca_free_annotation_predictions) <- paste("hlca", colnames(hlca_free_annotation_predictions), sep = "_")
cohort <- AddMetaData(cohort, metadata = hlca_free_annotation_predictions)
pdf(file="plots/hlca_umap_raw.pdf", height=6, width=10)
DimPlot(cohort, reduction = "umap", group.by = "hlca_predicted.id", label = TRUE,repel = TRUE, label.size=3) +NoLegend() 
dev.off()

#annotate cells with low prediction scores
hlcalabels=data.frame(Labels=cohort$hlca_predicted.id, Score=cohort$hlca_prediction.score.max)
hlcalabels$Labels[hlcalabels$Score<0.5]="Low_Score"
ourmeta=data.frame(hlcalabels=hlcalabels$Labels)
rownames(ourmeta)=rownames(hlcalabels)
cohort=AddMetaData(cohort,metadata = ourmeta)

pdf(file="plots/hlca_umap_all_lowscore.pdf", height=6, width=10)
DimPlot(cohort, reduction = "umap", group.by = "hlcalabels", label=T, pt.size = 0.0001, repel=T)
dev.off()

#make confusion matrix of scores
hlca=data.frame(cohort$merged, cohort$hlcalabels)
colnames(hlca)=c("Cluster","Label")
rownames(hlca)=NULL

hlcadata=data.frame(unclass(table(hlca)))
hlcamatrix=as.matrix(hlcadata)
rows=rowSums(hlcamatrix)
hlcapercent=hlcamatrix/rows
hlcapercent=t(hlcapercent)

pdf(file="plots/HLCA_confusionmatrix.pdf",height=8, width=10)
pheatmap(hlcapercent,cluster_rows = F, cluster_cols = F,
         color=viridis(100,option="viridis"),
         border_color =NA,
         angle_col = 45,cellwidth = 20, cellheight = 20,
         show_rownames = T, show_colnames = T)
dev.off()

##max in group
max=as.data.frame(hlcapercent)
max$CellType=rownames(hlcapercent)

clusternames=as.vector(unique(cohort$merged))
for(cluster in clusternames){
  individual=data.frame(ID=max[,cluster], CellType=max$CellType)
  sorted=individual[order(individual$ID, decreasing=T),]
  highest=sorted[1:2,]
  print(cluster)
  print(highest)
}

clustervect=as.vector(cohort$merged)
cohort$HLCATop=clustervect
cohort$HLCATop[cohort$HLCATop==1]="Monocyte"
cohort$HLCATop[cohort$HLCATop=="2_7"]="Monocyte"
cohort$HLCATop[cohort$HLCATop=="3_6"]="Monocyte"
cohort$HLCATop[cohort$HLCATop==4]="Monocyte"
cohort$HLCATop[cohort$HLCATop==5]="Monocyte"
cohort$HLCATop[cohort$HLCATop==8]="Differentiating Basal"
cohort$HLCATop[cohort$HLCATop==9]="Dendritic"
cohort$HLCATop[cohort$HLCATop==10]="Ciliated"
cohort$HLCATop[cohort$HLCATop==11]="Dendritic"

pdf(file="plots/hlca_relabeled_umaps.pdf", height=6, width=8)
DimPlot(cohort, reduction = "umap", group.by = "HLCATop", pt.size = 0.0001)
dev.off()

#saveRDS(cohort, file="Final_SeuratObjects/cohort_adams_hlca_labeltrans.RDS")


#####################################################################

#Adams label transfer

adams <-readRDS("/home/u30/hwelfley/adams_paper/adams_healthy.RDS")
adams <- NormalizeData(adams)
adams = FindVariableFeatures(adams, selection.method = "vst", nfeatures = 2000)
adams <- ScaleData(adams, verbose=FALSE)
adams <- RunPCA(adams, npcs = 30, verbose = FALSE)


# adams label transfer
adams.anchors <- FindTransferAnchors(reference = adams, query = cohort, dims = 1:30)
adams_predictions <- TransferData(anchorset = adams.anchors, refdata = adams$Manuscript_Identity, dims = 1:30)
colnames(adams_predictions) <- paste("adams", colnames(adams_predictions), sep = "_")
cohort <- AddMetaData(cohort, metadata = adams_predictions)
pdf(file="plots/adams_umap_raw.pdf", height=6, width=10)
DimPlot(cohort, reduction = "umap", group.by = "adams_predicted.id", label = TRUE,repel = TRUE, label.size=3) +NoLegend() 
dev.off()


adamlabels=data.frame(AdamsLabels=cohort$adams_predicted.id, Score=cohort$adams_prediction.score.max)
adamlabels$AdamsLabels[adamlabels$Score<0.5]="Low_Score"
ameta=data.frame(adamlabels=adamlabels$AdamsLabels)
rownames(ameta)=rownames(adamlabels)
cohort=AddMetaData(cohort,metadata = ameta)

pdf(file="plots/cohort_umap_all_lowscore.pdf", height=6, width=10)
DimPlot(cohort, reduction = "umap", group.by = "adamlabels", label=T, pt.size = 0.0001, repel=T)
dev.off()

#make confusion matrix of scores
adam=data.frame(cohort$merged, cohort$adamlabels)
colnames(adam)=c("Cluster","Label")
rownames(adam)=NULL

adamdata=data.frame(unclass(table(adam)))
adammatrix=as.matrix(adamdata)
arows=rowSums(adammatrix)
apercent=adammatrix/arows
apercent=t(apercent)

pdf(file="plots/adams_confusionmatrix.pdf",height=8, width=10)
pheatmap(apercent,cluster_rows = F, cluster_cols = F,
         color=viridis(100,option="magma"),
         border_color =NA,
         angle_col = 45,cellwidth = 20, cellheight = 20,
         show_rownames = T, show_colnames = T)
dev.off()

##max in group
amax=as.data.frame(apercent)
amax$CellType=rownames(apercent)

clusternames=as.vector(unique(cohort$merged))
for(cluster in clusternames){
  individual=data.frame(ID=amax[,cluster], CellType=amax$CellType)
  sorted=individual[order(individual$ID, decreasing=T),]
  highest=sorted[1:2,]
  print(cluster)
  print(highest)
}

aclustervect=as.vector(cohort$merged)
cohort$AdamsTop=aclustervect
cohort$AdamsTop[cohort$AdamsTop==1]="Monocyte"
cohort$AdamsTop[cohort$AdamsTop=="2_7"]="Macrophage"
cohort$AdamsTop[cohort$AdamsTop=="3_6"]="Monocyte"
cohort$AdamsTop[cohort$AdamsTop==4]="Macrophage"
cohort$AdamsTop[cohort$AdamsTop==5]="Macrophage"
cohort$AdamsTop[cohort$AdamsTop==8]="Goblet"
cohort$AdamsTop[cohort$AdamsTop==9]="Macrophage"
cohort$AdamsTop[cohort$AdamsTop==10]="Ciliated"
cohort$AdamsTop[cohort$AdamsTop==11]="Dendritic"

pdf(file="plots/adams_relabeled_umaps.pdf", height=6, width=10)
DimPlot(cohort, reduction = "umap", group.by = "AdamsTop", pt.size = 0.0001)
dev.off()

saveRDS(cohort,"cohort_adams_hlca.RDS")

###########################
#SingleR
#########################

library(SingleR)
library(celldex)
library(Seurat)
library(SeuratWrappers)

new.data=as.SingleCellExperiment(cohort)
saveRDS(new.data, file="cohort_sce.RDS")

ref.data <- HumanPrimaryCellAtlasData(ensembl=F)

predictions <- SingleR(test=new.data, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

saveRDS(predictions, file="cohort_singleR_predictions.RDS")

singlelabel=data.frame("SingleRLabels"=predictions$labels)
rownames(singlelabel)=predictions@rownames

cohort=AddMetaData(cohort, metadata = singlelabel)

pdf(file="plots/SingleR_raw_umap.pdf", height=6, width=8)
DimPlot(cohort, reduction="umap", group.by = "SingleRLabels", label=T, pt.size = 0.01)
dev.off()

labelinfo=data.frame(unclass(table(cohort$SingleRLabels,cohort$merged)))

cohort$TopSingleR=as.character(cohort$merged)
cohort$TopSingleR[cohort$TopSingleR==1]="Neutrophil"
cohort$TopSingleR[cohort$TopSingleR=="2_7"]="Monocyte"
cohort$TopSingleR[cohort$TopSingleR=="3_6"]="Monocyte"
cohort$TopSingleR[cohort$TopSingleR==4]="Monocyte"
cohort$TopSingleR[cohort$TopSingleR==5]="Monocyte"
cohort$TopSingleR[cohort$TopSingleR==8]="Epithelial"
cohort$TopSingleR[cohort$TopSingleR==9]="Monocyte"
cohort$TopSingleR[cohort$TopSingleR==10]="Epithelial"
cohort$TopSingleR[cohort$TopSingleR==11]="Monocyte"

pdf(file="plots/SingleR_top.pdf", height=6, width=8)
DimPlot(cohort, reduction="umap", group.by = "TopSingleR", pt.size = 0.01)
dev.off()

saveRDS(cohort, file="cohort_labeltrans_singler.RDS")

singleheat=data.frame(cohort$merged, cohort$SingleRLabels)
colnames(singleheat)=c("Cluster","Label")
rownames(singleheat)=NULL

singledata=data.frame(unclass(table(singleheat)))
singlematrix=as.matrix(singledata)
rrows=rowSums(singlematrix)
rpercent=singlematrix/rrows
rpercent=t(rpercent)

pdf(file="plots/SingleR_confusionmatrix.pdf",height=8, width=10)
pheatmap(rpercent,cluster_rows = F, cluster_cols = F,
         color=viridis(100,option="cividis"),
         border_color =NA,
         angle_col = 45,cellwidth = 20, cellheight = 20,
         show_rownames = T, show_colnames = T)
dev.off()


######################

#confusion matrices 

library(Seurat)
library(pheatmap)
library(viridis)

setwd("/Users/hollywelfley/Desktop/New_Aspirate/")
cohort=readRDS("Objects/cohort_labeltrans_singler.RDS")


#make confusion matrix of scores
hlca=data.frame(cohort$merged, cohort$hlcalabels)
colnames(hlca)=c("Cluster","Label")
rownames(hlca)=NULL

hlcadata=data.frame(unclass(table(hlca)))
hlcamatrix=as.matrix(hlcadata)
rows=rowSums(hlcamatrix)
hlcapercent=hlcamatrix/rows
hlcapercent=t(hlcapercent)

adam=data.frame(cohort$merged, cohort$adamlabels)
colnames(adam)=c("Cluster","Label")
rownames(adam)=NULL

adamdata=data.frame(unclass(table(adam)))
adammatrix=as.matrix(adamdata)
arows=rowSums(adammatrix)
apercent=adammatrix/arows
apercent=t(apercent)


singleheat=data.frame(cohort$merged, cohort$SingleRLabels)
colnames(singleheat)=c("Cluster","Label")
rownames(singleheat)=NULL

singledata=data.frame(unclass(table(singleheat)))
singlematrix=as.matrix(singledata)
rrows=rowSums(singlematrix)
rpercent=singlematrix/rrows
rpercent=t(rpercent)

write.table(hlcapercent, file="Supplement/LabelTransfer/hlca_labeltrans_confusionmatrix_values.txt", col.names = T, row.names = T)
write.table(apercent, file="Supplement/LabelTransfer/adams_labeltrans_confusionmatrix_values.txt", col.names = T, row.names = T)
write.table(rpercent, file="Supplement/LabelTransfer/singleR_labeltrans_confusionmatrix_values.txt", col.names = T, row.names = T)


