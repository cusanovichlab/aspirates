# aspirates
Code to perform the analysis described in Welfley et al.

Count matrices were generated with CellRanger and doublets were detected with `scrublet.py` \
Seurat objects with SoupX correction was performed with `seruatobjects_soupx_scrublet.R`  \
Generation of UMAP, clustering correlation, and cell type markers were found with `UMAP_clustering_markers_labeling.R` \
Cell type annotation was performed in `labeltransfer_singleR.R` 

Cell type by subject and sex analysis `pseudobulk_pca_spearmancor.R` \
Integration of neonatal and adult samples was performed with `adult_neonate_integration.R` \
Myeloid trajectory analysis generated with `trajectory_pseudotime.R`
