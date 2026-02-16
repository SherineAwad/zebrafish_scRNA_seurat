#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(ggplot2)
})

# ---- args (ONLY TWO) ----
parser <- ArgumentParser()
parser$add_argument("--rds", required = TRUE)
parser$add_argument("--markers", required = TRUE)
args <- parser$parse_args()

cat("Reading RDS...\n")
obj <- readRDS(args$rds)

cat("Reading markers...\n")
genes <- readLines(args$markers)
genes <- genes[genes != ""]

# ---- force assay ----
DefaultAssay(obj) <- "RNA"

# ---- check UMAP ----
if (!"umap" %in% Reductions(obj)) {
  stop("ERROR: No UMAP reduction named 'umap'")
}

# ---- build gene symbol map ----
rn <- rownames(obj)
symbols <- sub(".*~~", "", rn)   # keep part after ~~

gene_map <- setNames(rn, symbols)
# names(gene_map) = symbols
# values = full Seurat rownames

project <- tools::file_path_sans_ext(basename(args$rds))

cat("Total genes in object:", length(rn), "\n")

# ---- plot ----
for (g in genes) {

  if (!g %in% names(gene_map)) {
    cat("SKIP (symbol not found):", g, "\n")
    next
  }

  feature_id <- gene_map[[g]]
  cat("Plotting:", g, "->", feature_id, "\n")

  p <- FeaturePlot(
    object = obj,
    features = feature_id,
    reduction = "umap"
  ) + ggtitle(g)

  ggsave(
    filename = paste0(project, "_", g, ".png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

