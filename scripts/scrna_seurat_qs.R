#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

ensure_packages <- function(packages, auto_install = FALSE) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0 && isTRUE(auto_install)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  loaded <- vapply(packages, function(pkg) {
    suppressPackageStartupMessages(require(pkg, character.only = TRUE))
  }, FUN.VALUE = logical(1))
  if (!all(loaded)) {
    stop("Failed to load packages: ", paste(packages[!loaded], collapse = ", "))
  }
}

save_plot_multiple <- function(plot_obj, file_stem, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(filename = paste0(file_stem, ".png"), plot = plot_obj, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filename = paste0(file_stem, ".pdf"), plot = plot_obj, width = width, height = height)
}

detect_marker_fc_col <- function(markers_df) {
  if ("avg_log2FC" %in% colnames(markers_df)) return("avg_log2FC")
  if ("avg_logFC" %in% colnames(markers_df)) return("avg_logFC")
  # Seurat v5 may use "avg_log2FC"
  stop("Cannot detect logFC column in markers table. Columns found: ", paste(colnames(markers_df), collapse = ", "))
}

top_markers_per_cluster <- function(markers_df, n = 10) {
  fc_col <- detect_marker_fc_col(markers_df)
  split_list <- split(markers_df, markers_df$cluster)
  top_list <- lapply(split_list, function(df) {
    df <- df[order(df[[fc_col]], decreasing = TRUE), , drop = FALSE]
    utils::head(df, n)
  })
  res <- do.call(rbind, top_list)
  res
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to input .qs file (Seurat object or counts)", metavar = "file"),
  make_option(c("-o", "--outdir"), type = "character", default = "./scrna_output", help = "Output directory [default %default]"),
  make_option(c("--species"), type = "character", default = "human", help = "Species: human|mouse|custom [default %default]"),
  make_option(c("--mt-pattern"), type = "character", default = NA, help = "Custom mitochondrial gene regex (overrides species default)"),
  make_option(c("--ribo-pattern"), type = "character", default = NA, help = "Custom ribosomal gene regex (overrides species default)"),
  make_option(c("--min-features"), type = "integer", default = 200, help = "Min nFeature_RNA to keep cell [default %default]"),
  make_option(c("--min-counts"), type = "integer", default = 500, help = "Min nCount_RNA to keep cell [default %default]"),
  make_option(c("--max-mito"), type = "double", default = 20, help = "Max percent.mt to keep cell [default %default]"),
  make_option(c("--max-ribo"), type = "double", default = 60, help = "Max percent.ribo to keep cell [default %default]"),
  make_option(c("--normalize"), type = "character", default = "sct", help = "Normalization: sct|lognormalize [default %default]"),
  make_option(c("--regress-mt"), action = "store_true", default = TRUE, help = "Regress out percent.mt in normalization [default %default]"),
  make_option(c("--dims"), type = "integer", default = 30, help = "Number of PCs/UMAP dims [default %default]"),
  make_option(c("--resolutions"), type = "character", default = "0.2,0.4,0.6,0.8,1.0", help = "Comma-separated clustering resolutions [default %default]"),
  make_option(c("--top-markers"), type = "integer", default = 10, help = "Top markers per cluster to export [default %default]"),
  make_option(c("--min-pct"), type = "double", default = 0.25, help = "FindAllMarkers min.pct [default %default]"),
  make_option(c("--logfc-threshold"), type = "double", default = 0.25, help = "FindAllMarkers logfc.threshold [default %default]"),
  make_option(c("--seed"), type = "integer", default = 777, help = "Random seed [default %default]"),
  make_option(c("--n-cores"), type = "integer", default = 1, help = "Parallel workers if future is available [default %default]"),
  make_option(c("--auto-install"), action = "store_true", default = FALSE, help = "Auto-install missing CRAN packages [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  print_help(OptionParser(option_list = option_list))
  stop("--input is required")
}

outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(outdir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
markers_dir <- file.path(outdir, "markers")
dir.create(markers_dir, recursive = TRUE, showWarnings = FALSE)
logs_dir <- file.path(outdir, "logs")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "pipeline.log")
sink(log_file, split = TRUE)
on.exit({
  sink()
}, add = TRUE)

message("scrna_seurat_qs - starting")
message("Input: ", opt$input)
message("Outdir: ", outdir)

set.seed(opt$seed)

base_pkgs <- c("qs", "ggplot2", "patchwork")
seurat_pkgs <- c("Seurat")
optional_pkgs <- c("sctransform", "future", "Matrix")

ensure_packages(c(base_pkgs, seurat_pkgs, optional_pkgs), auto_install = opt$`auto-install`)

if (requireNamespace("future", quietly = TRUE) && opt$`n-cores` > 1) {
  future::plan(strategy = future::multisession, workers = opt$`n-cores`)
}

# Determine gene patterns
species <- tolower(opt$species)
mt_pat <- opt$`mt-pattern`
ribo_pat <- opt$`ribo-pattern`
if (is.na(mt_pat) || is.null(mt_pat) || mt_pat == "") {
  mt_pat <- switch(species,
    human = "^MT-",
    mouse = "^mt-",
    custom = "^MT-",
    "^MT-"
  )
}
if (is.na(ribo_pat) || is.null(ribo_pat) || ribo_pat == "") {
  ribo_pat <- switch(species,
    human = "^RP[SL]",
    mouse = "^Rp[sl]",
    custom = "^RP[SL]",
    "^RP[SL]"
  )
}
message("Using patterns: mito=", mt_pat, ", ribo=", ribo_pat)

# Read input .qs
obj_raw <- qs::qread(opt$input)
obj_class <- class(obj_raw)[1]
message("Read object of class: ", obj_class)

obj <- NULL
if (inherits(obj_raw, "Seurat")) {
  obj <- obj_raw
} else if (inherits(obj_raw, "SingleCellExperiment")) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required to convert SingleCellExperiment")
  obj <- Seurat::as.Seurat(obj_raw)
} else if (inherits(obj_raw, "dgCMatrix") || inherits(obj_raw, "Matrix") || is.matrix(obj_raw)) {
  obj <- Seurat::CreateSeuratObject(counts = obj_raw)
} else if (is.list(obj_raw) && !is.null(obj_raw$counts)) {
  obj <- Seurat::CreateSeuratObject(counts = obj_raw$counts)
} else {
  stop("Unsupported object type from .qs: ", paste(class(obj_raw), collapse = ", "))
}

DefaultAssay(obj) <- if ("RNA" %in% names(Assays(obj))) "RNA" else names(Assays(obj))[1]

# Compute QC metrics
message("Computing QC metrics...")
obj$percent.mt <- Seurat::PercentageFeatureSet(obj, pattern = mt_pat)
obj$percent.ribo <- Seurat::PercentageFeatureSet(obj, pattern = ribo_pat)

# QC plots
vln <- Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.1)
save_plot_multiple(vln, file.path(plots_dir, "qc_vln"), width = 10, height = 8)

fs1 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
fs2 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
save_plot_multiple(fs1 | fs2, file.path(plots_dir, "qc_scatter"), width = 12, height = 6)

# Filtering
message("Filtering cells...")
cells_keep <- Seurat::WhichCells(
  obj,
  expression = nFeature_RNA >= opt$`min-features` & nCount_RNA >= opt$`min-counts` & percent.mt <= opt$`max-mito` & percent.ribo <= opt$`max-ribo`
)
obj <- subset(obj, cells = cells_keep)
message("Cells kept: ", length(cells_keep))

# Normalization
normalize_method <- tolower(opt$normalize)
message("Normalization method: ", normalize_method)
if (normalize_method == "sct") {
  sct_vars <- if (isTRUE(opt$`regress-mt`)) c("percent.mt") else NULL
  sct_args <- list(object = obj, vars.to.regress = sct_vars, verbose = FALSE)
  if (requireNamespace("glmGamPoi", quietly = TRUE)) {
    sct_args$method <- "glmGamPoi"
  }
  obj <- do.call(Seurat::SCTransform, sct_args)
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = seq_len(opt$dims), reduction = "pca", verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = seq_len(opt$dims), verbose = FALSE)
} else if (normalize_method == "lognormalize") {
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, vars.to.regress = if (isTRUE(opt$`regress-mt`)) c("percent.mt") else NULL, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = seq_len(opt$dims), reduction = "pca", verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = seq_len(opt$dims), verbose = FALSE)
} else {
  stop("Unknown normalization method: ", normalize_method)
}

# PCA/UMAP plots
pca_elbow <- Seurat::ElbowPlot(obj, ndims = max(50, opt$dims))
save_plot_multiple(pca_elbow, file.path(plots_dir, "pca_elbow"), width = 6, height = 5)
umap_by_counts <- Seurat::FeaturePlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), reduction = "umap", ncol = 3)
save_plot_multiple(umap_by_counts, file.path(plots_dir, "umap_qc_features"), width = 12, height = 4)

# Clustering over resolutions
res_vec <- as.numeric(strsplit(opt$resolutions, ",")[[1]])
res_vec <- res_vec[is.finite(res_vec)]
message("Clustering at resolutions: ", paste(res_vec, collapse = ", "))
all_cluster_cols <- character(0)
for (res in res_vec) {
  obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
  # Try to find the correct metadata column name for clusters at this resolution
  candidate_cols <- c(
    paste0("SCT_snn_res.", res),
    paste0("RNA_snn_res.", res),
    paste0("integrated_snn_res.", res),
    paste0("seurat_clusters")
  )
  found <- candidate_cols[candidate_cols %in% colnames(obj@meta.data)]
  if (length(found) == 0) {
    # Fallback: try formatted numeric (e.g., 0.2 -> 0.2, 0.4 -> 0.4)
    res_chr <- format(res, trim = TRUE, scientific = FALSE)
    candidate_cols2 <- c(
      paste0("SCT_snn_res.", res_chr),
      paste0("RNA_snn_res.", res_chr)
    )
    found <- candidate_cols2[candidate_cols2 %in% colnames(obj@meta.data)]
  }
  if (length(found) == 0) {
    warning("Could not detect cluster column for resolution ", res, "; using Idents")
    cluster_col <- "seurat_clusters"
  } else {
    cluster_col <- found[1]
  }
  all_cluster_cols <- c(all_cluster_cols, cluster_col)
  # UMAP plot by cluster column
  p <- Seurat::DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE) + ggplot2::ggtitle(paste0("UMAP - res=", res))
  file_stem <- file.path(plots_dir, paste0("umap_res_", res))
  save_plot_multiple(p, file_stem, width = 7, height = 6)
}

# Marker detection per resolution
for (idx in seq_along(res_vec)) {
  res <- res_vec[idx]
  cluster_col <- all_cluster_cols[idx]
  if (!cluster_col %in% colnames(obj@meta.data)) {
    next
  }
  message("Finding markers for resolution ", res, " (", cluster_col, ") ...")
  Seurat::Idents(obj) <- obj@meta.data[[cluster_col]]
  markers <- Seurat::FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = opt$`min-pct`,
    logfc.threshold = opt$`logfc-threshold`,
    verbose = FALSE
  )
  markers_file <- file.path(markers_dir, paste0("markers_res_", res, ".csv"))
  utils::write.csv(markers, markers_file, row.names = FALSE)
  # Top markers
  top_df <- top_markers_per_cluster(markers, n = opt$`top-markers`)
  top_file <- file.path(markers_dir, paste0("top", opt$`top-markers`, "_markers_res_", res, ".csv"))
  utils::write.csv(top_df, top_file, row.names = FALSE)
}

# Save processed object and metadata
qs_out <- file.path(outdir, "seurat_processed.qs")
rds_out <- file.path(outdir, "seurat_processed.rds")
meta_out <- file.path(outdir, "cell_metadata.csv")
message("Saving outputs...")
qs::qsave(obj, qs_out, preset = "high")
saveRDS(obj, rds_out)
utils::write.csv(obj@meta.data, meta_out)

session_file <- file.path(outdir, "sessionInfo.txt")
utils::capture.output(utils::sessionInfo(), file = session_file)
message("Done. Outputs saved to ", outdir)

