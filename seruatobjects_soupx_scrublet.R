#soupx
library(SoupX)
library(Seurat)

setwd("/xdisk/darrenc/hwelfley/redo/91bp/without_testrun/")

####soupx
#load in 10x outs (needs to have filtered and raw matrix as well as cluster info)
thirteen_soup=load10X("NEA_13/outs/")
fifteen_soup=load10X("NEA_15/outs/")
twenty_soup=load10X("NEA_20/outs/")
twentytwo_soup=load10X("NEA_22/outs/")
twentythree_soup=load10X("NEA_23/outs/")
twentysix_soup=load10X("NEA_26/outs/")
twentyeight_soup=load10X("NEA_28/outs/")
twentynine_soup=load10X("NEA_29/outs/")
thirtyone_soup=load10X("NEA_31/outs/")
thirtythree_soup=load10X("NEA_33/outs/")

#estimate rho
thirteen_soup=autoEstCont(thirteen_soup, tfidfMin = 0.7)
#Estimated global rho of 0.22
fifteen_soup=autoEstCont(fifteen_soup, tfidfMin = 0.7)
#Estimated global rho of 0.01
twenty_soup=autoEstCont(twenty_soup, tfidfMin = 0.7)
#Estimated global rho of 0.01
twentytwo_soup=autoEstCont(twentytwo_soup, tfidfMin = 0.7)
#Estimated global rho of 0.05
twentythree_soup=autoEstCont(twentythree_soup, tfidfMin = 0.7)
#Estimated global rho of 0.27
twentysix_soup=autoEstCont(twentysix_soup, tfidfMin = 0.7)
#Estimated global rho of 0.04
twentyeight_soup=autoEstCont(twentyeight_soup, tfidfMin = 0.7)
#Estimated global rho of 0.01
twentynine_soup=autoEstCont(twentynine_soup, tfidfMin = 0.7)
#Estimated global rho of 0.01
thirtyone_soup=autoEstCont(thirtyone_soup, tfidfMin = 0.7)
#Estimated global rho of 0.09
thirtythree_soup=autoEstCont(thirtythree_soup, tfidfMin = 0.7)
#Estimated global rho of 0.01

#adjust
thirteen_out = adjustCounts(thirteen_soup)
fifteen_out = adjustCounts(fifteen_soup)
twenty_out = adjustCounts(twenty_soup)
twentytwo_out = adjustCounts(twentytwo_soup)
twentythree_out = adjustCounts(twentythree_soup)
twentysix_out = adjustCounts(twentysix_soup)
twentyeight_out = adjustCounts(twentyeight_soup)
twentynine_out = adjustCounts(twentynine_soup)
thirtyone_out = adjustCounts(thirtyone_soup)
thirtythree_out = adjustCounts(thirtythree_soup)

################# Seurat

#makeobjects,min feat 200
thirteen <-CreateSeuratObject(counts=thirteen_out, min.features = 200)
fifteen <-CreateSeuratObject(counts=fifteen_out, min.features = 200)
twenty <-CreateSeuratObject(counts=twenty_out, min.features = 200)
twentytwo <-CreateSeuratObject(counts=twentytwo_out, min.features = 200)
twentythree <-CreateSeuratObject(counts=twentythree_out, min.features = 200)
twentysix <-CreateSeuratObject(counts=twentysix_out, min.features = 200)
twentyeight <-CreateSeuratObject(counts=twentyeight_out, min.features = 200)
twentynine <-CreateSeuratObject(counts=twentynine_out, min.features = 200)
thirtyone <-CreateSeuratObject(counts=thirtyone_out, min.features = 200)
thirtythree <-CreateSeuratObject(counts=thirtythree_out, min.features = 200)


#subset on singlets determined by scrublet
#see python script scrublet.py
singlet_thirteen <- readLines("NEA_13/outs/thirteen_singlet_barcodes.csv")
singlet_fifteen <- readLines("NEA_15/outs/fifteen_singlet_barcodes.csv")
singlet_twenty <- readLines("NEA_20/outs/twenty_singlet_barcodes.csv")
singlet_twentytwo <- readLines("NEA_22/outs/twentytwo_singlet_barcodes.csv")
singlet_twentythree <- readLines("NEA_23/outs/twentythree_singlet_barcodes.csv")
singlet_twentysix <- readLines("NEA_26/outs/twentysix_singlet_barcodes.csv")
singlet_twentyeight <- readLines("NEA_28/outs/twentyeight_singlet_barcodes.csv")
singlet_twentynine <- readLines("NEA_29/outs/twentynine_singlet_barcodes.csv")
singlet_thirtyone <- readLines("NEA_31/outs/thirtyone_singlet_barcodes.csv")
singlet_thirtythree <- readLines("NEA_33/outs/thirtythree_singlet_barcodes.csv")


#add -1
fixed_singlet_thirteen <- paste(singlet_thirteen,"1",sep="-")
fixed_singlet_fifteen <- paste(singlet_fifteen,"1",sep="-")
fixed_singlet_twenty <- paste(singlet_twenty,"1",sep="-")
fixed_singlet_twentytwo <- paste(singlet_twentytwo,"1",sep="-")
fixed_singlet_twentythree <- paste(singlet_twentythree,"1",sep="-")
fixed_singlet_twentysix <- paste(singlet_twentysix,"1",sep="-")
fixed_singlet_twentyeight <- paste(singlet_twentyeight,"1",sep="-")
fixed_singlet_twentynine <- paste(singlet_twentynine,"1",sep="-")
fixed_singlet_thirtyone <- paste(singlet_thirtyone,"1",sep="-")
fixed_singlet_thirtythree <- paste(singlet_thirtythree,"1",sep="-")

#subset on singlets
thirteen <- subset(thirteen, cells = fixed_singlet_thirteen)
fifteen<- subset(fifteen, cells = fixed_singlet_fifteen)
twenty<- subset(twenty, cells = fixed_singlet_twenty)
twentytwo<- subset(twentytwo, cells = fixed_singlet_twentytwo)
twentythree<- subset(twentythree, cells = fixed_singlet_twentythree)
twentysix <- subset(twentysix, cells = fixed_singlet_twentysix)
twentyeight <- subset(twentyeight, cells = fixed_singlet_twentyeight)
twentynine<- subset(twentynine, cells = fixed_singlet_twentynine)
thirtyone<- subset(thirtyone, cells = fixed_singlet_thirtyone)
thirtythree<- subset(thirtythree, cells = fixed_singlet_thirtythree)


#add percent mt and subset on mt filter
thirteen[["percent.mt"]] <- PercentageFeatureSet(thirteen, pattern = "^MT-")
fifteen[["percent.mt"]] <- PercentageFeatureSet(fifteen, pattern = "^MT-")
twenty[["percent.mt"]] <- PercentageFeatureSet(twenty, pattern = "^MT-")
twentytwo[["percent.mt"]] <- PercentageFeatureSet(twentytwo, pattern = "^MT-")
twentythree[["percent.mt"]] <- PercentageFeatureSet(twentythree, pattern = "^MT-")
twentysix[["percent.mt"]] <- PercentageFeatureSet(twentysix, pattern = "^MT-")
twentyeight[["percent.mt"]] <- PercentageFeatureSet(twentyeight, pattern = "^MT-")
twentynine[["percent.mt"]] <- PercentageFeatureSet(twentynine, pattern = "^MT-")
thirtyone[["percent.mt"]] <- PercentageFeatureSet(thirtyone, pattern = "^MT-")
thirtythree[["percent.mt"]] <- PercentageFeatureSet(thirtythree, pattern = "^MT-")


###subset on less than 20 percent mt reads
thirteen <- subset(thirteen, subset = percent.mt < 20)
fifteen<- subset(fifteen, subset = percent.mt < 20)
twenty <- subset(twenty, subset = percent.mt < 20)
twentytwo <- subset(twentytwo, subset = percent.mt < 20)
twentythree <- subset(twentythree, subset = percent.mt < 20)
twentysix<- subset(twentysix, subset = percent.mt < 20)
twentyeight<- subset(twentyeight, subset = percent.mt < 20)
twentynine<- subset(twentynine, subset = percent.mt < 20)
thirtyone<- subset(thirtyone, subset = percent.mt < 20)
thirtythree<- subset(thirtythree, subset = percent.mt < 20)


#add orig.ident
thirteen@meta.data$orig.ident="NEA_13"
fifteen$orig.ident="NEA_15"
twenty@meta.data$orig.ident="NEA_20"
twentytwo@meta.data$orig.ident="NEA_22"
twentythree@meta.data$orig.ident="NEA_23"
twentysix@meta.data$orig.ident="NEA_26"
twentyeight@meta.data$orig.ident="NEA_28"
twentynine$orig.ident="NEA_29"
thirtyone@meta.data$orig.ident="NEA_31"
thirtythree$orig.ident="NEA_33"

#add sex
thirteen$sex="M"
fifteen$sex="F"
twenty$sex="M"
twentytwo$sex="M"
twentythree$sex="F"
twentysix$sex="F"
twentyeight$sex="M"
twentynine$sex="F"
thirtyone$sex="M"
thirtythree$sex="M"

#add batch
thirteen$batch="B1"
fifteen$batch="B1"
twenty$batch="B1"
twentytwo$batch="B1"
twentythree$batch="B1"
twentysix$batch="B2"
twentyeight$batch="B2"
twentynine$batch="B2"
thirtyone$batch="B2"
thirtythree$batch="B2"

#saveobjects
saveRDS(thirteen, file="Final_SeuratObjects/NEA_13_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(fifteen, file="Final_SeuratObjects/NEA_15_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(twenty, file="Final_SeuratObjects/NEA_20_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(twentytwo, file="Final_SeuratObjects/NEA_22_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(twentythree, file="Final_SeuratObjects/NEA_23_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(twentysix, file="Final_SeuratObjects/NEA_26_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(twentyeight, file="Final_SeuratObjects/NEA_28_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(twentynine, file="Final_SeuratObjects/NEA_29_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(thirtyone, file="Final_SeuratObjects/NEA_31_seuratobj_soupx_scrublet200minfeatmt20.RDS")
saveRDS(thirtythree, file="Final_SeuratObjects/NEA_33_seuratobj_soupx_scrublet200minfeatmt20.RDS")

cohort=merge(thirteen, y=c(fifteen, twenty, twentytwo, twentythree, twentysix, twentyeight, twentynine, thirtyone, thirtythree),
      add.cell.ids=c("NEA_13","NEA_15","NEA_20","NEA_22","NEA_23","NEA_26","NEA_28","NEA_29","NEA_31","NEA_33"))

saveRDS(cohort, file="Final_SeuratObjects/cohort_onlymerged_seuratobj_soupx_scrublet200minfeatmt20.RDS")


##########################################################
