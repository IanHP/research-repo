
# 003 Simple Clustering
Using features selected in 002 to represent in 2D and cluster frames. Simple, as it is just skikit-learn.
Using specifically UMAP, PCA, and HDBSCAN

## Figures
This holds all figures output by below scripts

## Scripts
import_data.py: helper scripts to extract data from distance parquets
umap_graph: creates a UMAP dim reduction graph
pca_graph: creates a PCA dim reduction graph
hdbscan_cluster: will cluster together frames using hdbscan (unused)
centroid_extraction_tcl_creator: creates a tcl script that extracts centroid frames
    highlighted by hdbscan clustering
