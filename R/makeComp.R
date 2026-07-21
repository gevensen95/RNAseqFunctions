#' Run a differential expression comparison
#'
#' Takes a raw counts matrix and sample metadata, runs differential
#' expression with either DESeq2 or edgeR, and returns a results
#' data.frame shaped for direct use with [eRupt2()] and
#' [WebGestalt_Pipe2()]: `log2FoldChange` in column 2, gene IDs in column 7
#' (named "Gene"), and columns named "Gene", "padj", and "log2FoldChange".
#'
#' DESeq2 and edgeR are not hard dependencies of RNAseqBasic (to keep
#' installs light) — install whichever engine you plan to use via
#' `BiocManager::install("DESeq2")` or `BiocManager::install("edgeR")`.
#' Using a `contrast` with `method = "edgeR"` additionally requires the
#' `limma` package.
#'
#' @param counts Matrix or data.frame of raw (un-normalized) integer counts,
#' genes as rows, samples as columns. Row names must be gene identifiers.
#' @param metadata Data.frame of sample metadata. `rownames(metadata)` must
#' match `colnames(counts)` exactly, in the same order.
#' @param formula Design formula, e.g. `~condition`. If `NULL` (the
#' default), a `~condition` design is used provided `metadata` has a column
#' named "condition"; otherwise an error is raised asking for an explicit
#' formula.
#' @param contrast The comparison to extract.
#' - For `method = "DESeq2"`: a character vector of length 3, e.g.
#'   `c("condition", "treated", "control")`, passed to
#'   [DESeq2::results()]. If `NULL`, DESeq2's default contrast (the last
#'   level of the last formula variable vs. the reference level) is used.
#' - For `method = "edgeR"`: a single contrast string suitable for
#'   [limma::makeContrasts()], e.g. `"conditiontreated - conditioncontrol"`,
#'   evaluated against `model.matrix(formula, metadata)`. If `NULL`, the
#'   last coefficient of the design matrix is tested.
#' @param method Which engine to use: `"DESeq2"` (default) or `"edgeR"`.
#' @param alpha Significance cut-off passed through to DESeq2's
#' `results(alpha = )` for independent filtering. Ignored for edgeR.
#' @return A data.frame of DEG results with `log2FoldChange` as column 2
#' and gene IDs (as column "Gene") as column 7, ready to pass to
#' [eRupt2()] or [WebGestalt_Pipe2()].
#' @export
makeComp <- function(counts, metadata, formula = NULL, contrast = NULL,
                     method = c("DESeq2", "edgeR"), alpha = 0.05) {

  method <- match.arg(method)

  if (is.null(formula)) {
    if ("condition" %in% colnames(metadata)) {
      formula <- stats::as.formula("~condition")
    } else {
      stop("`formula` was not supplied and `metadata` has no 'condition' ",
           "column. Either pass a design formula (e.g. ~condition) or add ",
           "a 'condition' column to metadata.")
    }
  }
  formula <- stats::as.formula(formula)

  if (!identical(colnames(counts), rownames(metadata))) {
    stop("colnames(counts) must match rownames(metadata) exactly, ",
         "in the same order.")
  }

  if (method == "DESeq2") {
    result <- makeComp_DESeq2(counts, metadata, formula, contrast, alpha)
  } else {
    result <- makeComp_edgeR(counts, metadata, formula, contrast)
  }

  result
}

makeComp_DESeq2 <- function(counts, metadata, formula, contrast, alpha) {

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for method = \"DESeq2\". ",
         "Install it with BiocManager::install(\"DESeq2\").")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = metadata,
    design = formula
  )
  dds <- DESeq2::DESeq(dds)

  res <- if (is.null(contrast)) {
    DESeq2::results(dds, alpha = alpha)
  } else {
    DESeq2::results(dds, contrast = contrast, alpha = alpha)
  }

  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  rownames(res_df) <- NULL
  res_df
}

makeComp_edgeR <- function(counts, metadata, formula, contrast) {

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for method = \"edgeR\". ",
         "Install it with BiocManager::install(\"edgeR\").")
  }

  design_matrix <- stats::model.matrix(formula, data = metadata)

  y <- edgeR::DGEList(counts = round(as.matrix(counts)))
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design_matrix)
  fit <- edgeR::glmQLFit(y, design_matrix)

  if (is.null(contrast)) {
    qlf <- edgeR::glmQLFTest(fit, coef = ncol(design_matrix))
  } else {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Package 'limma' is required to use `contrast` with ",
           "method = \"edgeR\". Install it with ",
           "BiocManager::install(\"limma\").")
    }
    contrast_vec <- limma::makeContrasts(contrasts = contrast, levels = design_matrix)
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
  }

  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table

  stat_col <- if ("F" %in% colnames(tt)) tt$F else tt$LR

  res_df <- data.frame(
    baseMean = 2 ^ tt$logCPM,
    log2FoldChange = tt$logFC,
    lfcSE = NA_real_,
    stat = stat_col,
    pvalue = tt$PValue,
    padj = tt$FDR,
    Gene = rownames(tt),
    row.names = NULL
  )
  res_df
}
