# Internal helpers shared by PCAplot() and DEGheatmap() for turning raw
# counts into a normalized (log-scale) expression matrix. Not exported.

normalize_counts_deseq2 <- function(counts, metadata) {
  if (!requireNamespace("DESeq2", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Packages 'DESeq2' and 'SummarizedExperiment' are required for ",
         "method = \"DESeq2\". Install them with ",
         "BiocManager::install(c(\"DESeq2\", \"SummarizedExperiment\")).")
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = metadata,
    design = ~1
  )

  # vst() fits a parametric mean-dispersion trend using a subsample of
  # `nsub` genes (1000 by default) and errors out below that. Fall back to
  # the exact (slower, but always-applicable) transformation for smaller
  # gene sets, e.g. targeted panels or already-filtered results.
  if (nrow(dds) < 1000) {
    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  } else {
    vsd <- DESeq2::vst(dds, blind = TRUE)
  }

  SummarizedExperiment::assay(vsd)
}

normalize_counts_edger <- function(counts) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for method = \"edgeR\". ",
         "Install it with BiocManager::install(\"edgeR\").")
  }
  y <- edgeR::DGEList(counts = round(as.matrix(counts)))
  y <- edgeR::calcNormFactors(y)
  edgeR::cpm(y, log = TRUE)
}

normalize_counts <- function(counts, metadata, method) {
  switch(
    method,
    DESeq2 = normalize_counts_deseq2(counts, metadata),
    edgeR = normalize_counts_edger(counts)
  )
}
