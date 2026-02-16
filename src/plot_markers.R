#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(ggplot2)
})

# ---------------- args ----------------
parser <- ArgumentParser()
parser$add_argument("--rds", required = TRUE)
parser$add_argument("--markers", required = TRUE)
args <- parser$parse_args()

# ---------------- load ----------------
obj <- readRDS(args$rds)
genes <- readLines(args$markers)
genes <- genes[genes != ""]

DefaultAssay(obj) <- "RNA"

if (!"umap" %in% Reductions(obj)) {
  stop("ERROR: UMAP not found")
}

project <- tools::file_path_sans_ext(basename(args$rds))
dir.create("figures", showWarnings = FALSE)

# ---------------- ANNOTATED UMAP ----------------
# Use the actual column in your object
annotation_col <- "new_celltype"

p_umap <- DimPlot(
  object = obj,
  reduction = "umap",
  group.by = annotation_col,
  label = TRUE,
  repel = TRUE
)

ggsave(
  filename = file.path("figures", paste0(project, "_UMAP_", annotation_col, ".png")),
  plot = p_umap,
  width = 7,
  height = 6,
  dpi = 300
)

# ---------------- gene symbol mapping ----------------
rn <- rownames(obj)
symbols <- sub(".*~~", "", rn)
gene_map <- setNames(rn, symbols)

# ---------------- feature UMAPs ----------------
for (g in genes) {

  if (!g %in% names(gene_map)) {
    cat("SKIP:", g, "\n")
    next
  }

  p <- FeaturePlot(
    object = obj,
    features = gene_map[[g]],
    reduction = "umap"
  ) + ggtitle(g)

  ggsave(
    filename = file.path("figures", paste0(project, "_", g, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

