library(Seurat)

setwd("/Users/hollywelfley/Desktop/New_Aspirate/")

aspirate=readRDS("Objects/cohort_labeltrans_singler.RDS")

#split by cell type
basal=subset(aspirate, FinalLabel=="Basal")
ciliated=subset(aspirate, FinalLabel=="Ciliated")
cycling=subset(aspirate, FinalLabel=="Cycling")
ione=subset(aspirate, FinalLabel=="Intermediate I (VDAC1)")
itwo=subset(aspirate, FinalLabel=="Intermediate II (MIR222HG)")
ithree=subset(aspirate, FinalLabel=="Intermediate III (BEST1)")
monocyte=subset(aspirate, FinalLabel=="Monocyte")
macrophage=subset(aspirate, FinalLabel=="Macrophage")
neutrophil=subset(aspirate, FinalLabel=="Neutrophil")

#pseudobulk Seurat objects by subject
basalbulk=AverageExpression(basal, return.seurat = T, group.by = "orig.ident")
basalbulk$orig.ident=rownames(basalbulk@meta.data)
basalbulk$FinalLabel="Basal"

ciliatedbulk=AverageExpression(ciliated, return.seurat = T, group.by = "orig.ident")
ciliatedbulk$orig.ident=rownames(ciliatedbulk@meta.data)
ciliatedbulk$FinalLabel="Ciliated"

cyclingbulk=AverageExpression(cycling, return.seurat = T, group.by = "orig.ident")
cyclingbulk$orig.ident=rownames(cyclingbulk@meta.data)
cyclingbulk$FinalLabel="Cycling"

ionebulk=AverageExpression(ione, return.seurat = T, group.by = "orig.ident")
ionebulk$orig.ident=rownames(ionebulk@meta.data)
ionebulk$FinalLabel="Intermediate I (VDAC1)"

itwobulk=AverageExpression(itwo, return.seurat = T, group.by = "orig.ident")
itwobulk$orig.ident=rownames(itwobulk@meta.data)
itwobulk$FinalLabel="Intermediate II (MIR222HG)"

ithreebulk=AverageExpression(ithree, return.seurat = T, group.by = "orig.ident")
ithreebulk$orig.ident=rownames(ithreebulk@meta.data)
ithreebulk$FinalLabel="Intermediate III (BEST1)"

monocytebulk=AverageExpression(monocyte, return.seurat = T, group.by = "orig.ident")
monocytebulk$orig.ident=rownames(monocytebulk@meta.data)
monocytebulk$FinalLabel="Monocyte"

macrophagebulk=AverageExpression(macrophage, return.seurat = T, group.by = "orig.ident")
macrophagebulk$orig.ident=rownames(macrophagebulk@meta.data)
macrophagebulk$FinalLabel="Macrophage"

neutrophilbulk=AverageExpression(neutrophil, return.seurat = T, group.by = "orig.ident")
neutrophilbulk$orig.ident=rownames(neutrophilbulk@meta.data)
neutrophilbulk$FinalLabel="Neutrophil"

bulk=merge(x=basalbulk, y=c(ciliatedbulk, cyclingbulk, ionebulk, itwobulk, ithreebulk, monocytebulk, macrophagebulk, neutrophilbulk))
bulk$Sex=bulk$orig.ident
bulk$Sex[bulk$Sex %in% c("NEA_15","NEA_23","NEA_26","NEA_29")]="F"
bulk$Sex[bulk$Sex!="F"]="M"

#pca
bulk <-FindVariableFeatures(bulk, selection.method = "vst", nfeatures = 3000)
bulk <- ScaleData(bulk, verbose = FALSE)
bulk <- RunPCA(bulk, pc.genes = bulk@var.genes, npcs = 20, verbose = FALSE)

p1=DimPlot(bulk, reduction = "pca", group.by = "FinalLabel", pt.size = 1)
p2=DimPlot(bulk, reduction = "pca", group.by = "orig.ident", pt.size = 1)
p3=DimPlot(bulk, reduction = "pca", group.by = "Sex", pt.size = 1)

pdf(file="pseudobulk_pca.pdf", height=6, width=20)
p1+p2+p3
dev.off()

bulk$Sex=factor(bulk$Sex, levels=c("M","F"))
thisone=c("#008ECE","#59C7EB","#8E2043","#E0607E","#077187","#54BFB7","#3E5496","#8290BB","#ECA0B2")
pdf(file="pseudobulk_pca12_celltype_sexshape_bigger.pdf", height=6, width=8)
DimPlot(bulk, reduction = "pca", group.by = "FinalLabel", pt.size = 4, shape.by = "Sex", cols=thisone)
dev.off()

saveRDS(bulk, file="Objects/pseudobulk_by_sub_cell.RDS")

# Get the total variance:
mat <- GetAssayData(bulk, assay = "RNA", slot = "scale.data")
pca <- bulk[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

############################
#Spearman correlation
############################

#pseudobulk matrices by subject
basalbulkmatrix=AverageExpression(basal, return.seurat = F, group.by = "orig.ident")
basalbulkmatrix_na = basalbulkmatrix[[1]]
basalbulkmatrix_na[basalbulkmatrix_na == 0] = NA
colnames(basalbulkmatrix_na)=paste("Basal",colnames(basalbulkmatrix_na), sep="_")

ciliatedbulkmatrix=AverageExpression(ciliated, return.seurat = F, group.by = "orig.ident")
ciliatedbulkmatrix_na = ciliatedbulkmatrix[[1]]
ciliatedbulkmatrix_na[ciliatedbulkmatrix_na == 0] = NA
colnames(ciliatedbulkmatrix_na)=paste("Ciliated",colnames(ciliatedbulkmatrix_na), sep="_")

cyclingbulkmatrix=AverageExpression(cycling, return.seurat = F, group.by = "orig.ident")
cyclingbulkmatrix_na = cyclingbulkmatrix[[1]]
cyclingbulkmatrix_na[cyclingbulkmatrix_na == 0] = NA
colnames(cyclingbulkmatrix_na)=paste("Cycling",colnames(cyclingbulkmatrix_na), sep="_")

ionebulkmatrix=AverageExpression(ione, return.seurat = F, group.by = "orig.ident")
ionebulkmatrix_na = ionebulkmatrix[[1]]
ionebulkmatrix_na[ionebulkmatrix_na == 0] = NA
colnames(ionebulkmatrix_na)=paste("InterI",colnames(ionebulkmatrix_na), sep="_")

itwobulkmatrix=AverageExpression(itwo, return.seurat = F, group.by = "orig.ident")
itwobulkmatrix_na = itwobulkmatrix[[1]]
itwobulkmatrix_na[itwobulkmatrix_na == 0] = NA
colnames(itwobulkmatrix_na)=paste("InterII",colnames(itwobulkmatrix_na), sep="_")

ithreebulkmatrix=AverageExpression(ithree, return.seurat = F, group.by = "orig.ident")
ithreebulkmatrix_na = ithreebulkmatrix[[1]]
ithreebulkmatrix_na[ithreebulkmatrix_na == 0] = NA
colnames(ithreebulkmatrix_na)=paste("InterIII",colnames(ithreebulkmatrix_na), sep="_")

monocytebulkmatrix=AverageExpression(monocyte, return.seurat = F, group.by = "orig.ident")
monocytebulkmatrix_na = monocytebulkmatrix[[1]]
monocytebulkmatrix_na[monocytebulkmatrix_na == 0] = NA
colnames(monocytebulkmatrix_na)=paste("Monocyte",colnames(monocytebulkmatrix_na), sep="_")

macrophagebulkmatrix=AverageExpression(macrophage, return.seurat = F, group.by = "orig.ident")
macrophagebulkmatrix_na = macrophagebulkmatrix[[1]]
macrophagebulkmatrix_na[macrophagebulkmatrix_na == 0] = NA
colnames(macrophagebulkmatrix_na)=paste("Macrophage",colnames(macrophagebulkmatrix_na), sep="_")

neutrophilbulkmatrix=AverageExpression(neutrophil, return.seurat = F, group.by = "orig.ident")
neutrophilbulkmatrix_na = neutrophilbulkmatrix[[1]]
neutrophilbulkmatrix_na[neutrophilbulkmatrix_na == 0] = NA
colnames(neutrophilbulkmatrix_na)=paste("Neutrophil",colnames(neutrophilbulkmatrix_na), sep="_")

####

bulkmatrix=cbind(basalbulkmatrix_na,ciliatedbulkmatrix_na, cyclingbulkmatrix_na, ionebulkmatrix_na, itwobulkmatrix_na, ithreebulkmatrix_na, monocytebulkmatrix_na, macrophagebulkmatrix_na, neutrophilbulkmatrix_na)
corbulkmatrix=cor(bulkmatrix,method="spearman",use="pairwise.complete.obs")
write.table(corbulkmatrix, file="bulk_allcelltypes_subject_spearman_table.txt", col.names = T, row.names = F, quote = F)

## generate figure

subcelldata=data.frame("Subject_Cell"=colnames(corbulkmatrix))
subcelldata$CellType=sub("_NEA_.*","",subcelldata$Subject_Cell)
subcelldata$Subject=sub(".*_", "", subcelldata$Subject_Cell)
subcelldata$Sex=subcelldata$Subject
subcelldata$Sex[subcelldata$Sex %in% c("15","23","26","29")]="F"
subcelldata$Sex[subcelldata$Sex %in% c("13","20","22","28","31","33")]="M"

sexcelldata=data.frame("CellType"=subcelldata$CellType, "Sex"=subcelldata$Sex)
rownames(sexcelldata)=subcelldata$Subject_Cell

colorsexlist=list("Sex"=c("M"="#9AA0A7","F"="#FEA090"),
                  "CellType"=c("Basal"="#008ECE", "Ciliated"="#59C7EB","Cycling"="#8E2043",
                               "InterI"="#E0607E","InterII"="#077187", "InterIII"="#54BFB7",
                               "Monocyte"="#3E5496","Macrophage"="#8290BB", "Neutrophil"="#ECA0B2"))


pdf(file="Recent_Figures/bulk_allcelltypes_subject_spearman_sexcolorlab_clustered.pdf", height=8, width=10)
pheatmap(corbulkmatrix, col = brewer.pal(9,"Blues"), 
         annotation_col = sexcelldata, annotation_colors = colorsexlist,
         show_rownames = F, show_colnames = F, border_color = NA)
dev.off()

