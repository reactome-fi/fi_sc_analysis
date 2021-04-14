# ------------------------------------------
# Read args
# ------------------------------------------
args = commandArgs(trailingOnly=T)

if(length(args) < 2) {
  message("Invalid number of passed arguments.")
}

umi.path <- args[1] 
species <- args[2]
thresholds <- args[3]
n <- args[4]
k <- args[5]
plotting <- args[6]
brewer.name <- args[7]
# ------------------------------------------
# seed 
# ------------------------------------------
set.seed(1234)

# ------------------------------------------
# install UniPath via GitHub 
# https://reggenlab.github.io/UniPathWeb/
# library(devtools)
# install_github("reggenlab/UniPath")
# ------------------------------------------
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")
# ------------------------------------------
# install dependencies 
# ------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("netbiov")
# BiocManager::install("GenomicRanges")

p_load("vegan")
p_load("FNN")
p_load("igraph")
p_load("preprocessCore")
p_load("GenomicRanges")

p_load("netbiov")
p_load("GenomicRanges")

p_load("RColorBrewer")
# --------------------------------------------
# https://anndata.dynverse.org/index.html
p_load("anndata")
p_load("reticulate")
reticulate::use_python("/opt/anaconda3/bin/python", required = TRUE) # python path 
reticulate::py_config()

# --------------------------------------------
# load and install unipath from GitHub
# ------------------------------------------
p_load_gh("reggenlab/UniPath")

# --------------------------------------------
# user UMI 
# --------------------------------------------
ad <- read_h5ad("E17_adult_anndata.h5ad")
umi_expression <- t(as.data.frame(as.matrix(ad$X)))
species <- "mouse"
threshold <- 3 # genesets having number of genes greater than the threshold value provided
n <- 4 # number of clusters corresponding to type of cells 
k <- 5 # top k nearest neighbor computation
plotting <- T
color.brewer.name <- "Set2"

# --------------------------------------------
# mouse/human 
# load null model data matrix 
# load symbols/markers
# --------------------------------------------
if (species == "mouse"){
  
  data("mouse_null_model")
  data("c5.bp.v6.0.symbols")
  # ---------------------------------------
  # browns method to combine p-values of null model data matrix (pre-annotated)
  # --------------------------------------
  # message("Combining p-values...")
  Pval <- binorm(mouse_null_data) 
  combp_ref <- combine(c5.bp.v6.0.symbols, mouse_null_data, rownames(mouse_null_data), Pval, thr=threshold)
  
  # ---------------------------------------
  # User-defined expression 
  # --------------------------------------
  Pval1 <- binorm(umi_expression)
  combp <- combine(c5.bp.v6.0.symbols, umi_expression ,rownames(umi_expression), Pval1, thr=threshold)
  
} else if (species == "human") {
  
  data("human_null_model")
  data("human_markers")
  # ---------------------------------------
  # browns method to combine p-values of null model data matrix (pre-annotated)
  # --------------------------------------
  # message("Combining p-values...")
  Pval <- binorm(human_null_data)
  combp_ref <- combine(human_markers, human_null_data, rownames(human_null_data), Pval, thr=threshold)
  # ---------------------------------------
  # User-defined expression 
  # --------------------------------------
  Pval1 <- binorm(umi_expression)
  combp <- combine(human_markers, umi_expression ,rownames(umi_expression), Pval1, thr=threshold)
  
} else {
  
  message("Provide a species of interest.")
  
}

# ---------------------------------------
# The adjusted p-value matrix (scores$adjpvalog) is referred to as pathway scores.
# --------------------------------------
scores <- adjust(combp, combp_ref)
# save(scores, file = "scores.RData")
load("/Users/sanati/Documents/reactome_embeddings/scores.RData")

# ---------------------------------------
# Pseudo temporal ordering
# TODO: save/return results
# ---------------------------------------
distclust <- dist_clust(scores$adjpvalog, n=n)
dist <- distclust$distance
clusters <- distclust$clusters # cell clusters 
index <- index(scores$adjpvalog, k=k)
KNN <- KNN(scores$adjpvalog, index, clusters)
node_class <- class1(clusters, KNN)
distance <- distance(dist, node_class, clusters)
corr_mst <- minimum_spanning_tree(distance) # igraph object


# ---------------------------------------
# plotting
# ---------------------------------------
if (plotting == T){
  vertex_color <- brewer.pal(n = n, name = color.brewer.name)
  # TODO:  fetch cell_labels from anndata
  cell_labels <- data.frame(c(rep("E18.5",82), rep("E14.5",44), rep("Adult",46), rep("E16.5",23))) 
  # mst.plot.mod(corr_mst, vertex.color = vertex_color, mst.edge.col="black", 
  #                       bg="white", layout.function="layout.kamada.kawai")
  # Note: bug fix but edges don't draw (ok to move on)
  UniPath::mst.plot.mod(corr_mst, vertex.color = vertex_color[as.factor(cell_labels[,1])], mst.edge.col="black", bg="white", layout.function="layout.kamada.kawai", v.size = 3, e.size=0.005, mst.e.size = 0.005)
  legend("top", legend = sort(unique(cell_labels[,1])), col = vertex_color,pch=20, box.lty=0, cex=0.6, pt.cex=1.5, horiz=T)
  # TODO: save/return plot 
}

