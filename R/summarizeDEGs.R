#' Summarize DEG counts by significance
#'
#' Tallies up- and down-regulated genes from one or more DEG results
#' data.frames (e.g. from [makeComp()]), for a quick summary table across
#' comparisons.
#'
#' @param results A single results data.frame, or a named list of results
#' data.frames (e.g. from running [makeComp()] once per contrast).
#' @param alpha padj cut-off for significance
#' @param FCsig absolute log2FC cut-off for significance
#' @return A data.frame with one row per comparison and columns
#' "comparison", "up", "down", "not_significant", and "total".
#' @export
summarizeDEGs <- function(results, alpha = 0.05, FCsig = 0) {

  if (is.data.frame(results) || !is.list(results)) {
    results <- list(comparison = results)
  }

  if (is.null(names(results)) || any(names(results) == "")) {
    names(results) <- paste0("comparison_", seq_along(results))
  }

  rows <- lapply(names(results), function(nm) {
    df <- stats::na.omit(results[[nm]][, c("padj", "log2FoldChange")])
    up <- sum(df$padj < alpha & df$log2FoldChange > FCsig)
    down <- sum(df$padj < alpha & df$log2FoldChange < (FCsig * -1))
    data.frame(
      comparison = nm,
      up = up,
      down = down,
      not_significant = nrow(df) - up - down,
      total = nrow(df),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
