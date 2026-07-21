#' Sample-to-sample distance heatmap
#'
#' Normalizes raw counts and plots a heatmap of Euclidean distances between
#' samples (hierarchically clustered), for QC of sample similarity and
#' outlier detection. A complement to [PCAplot()].
#'
#' @param counts Matrix or data.frame of raw counts, genes as rows, samples
#' as columns (same orientation as [makeComp()]).
#' @param metadata Data.frame of sample metadata. `rownames(metadata)` must
#' match `colnames(counts)` exactly, in the same order.
#' @param method Which engine to use for normalization: `"DESeq2"`
#' (default, via `vst()`) or `"edgeR"` (via `cpm(log = TRUE)`).
#' @return A ggplot2 heatmap object.
#' @export
sampleDistHeatmap <- function(counts, metadata, method = c("DESeq2", "edgeR")) {

  method <- match.arg(method)

  if (!identical(colnames(counts), rownames(metadata))) {
    stop("colnames(counts) must match rownames(metadata) exactly, ",
         "in the same order.")
  }

  norm_mat <- normalize_counts(counts, metadata, method)

  sample_dist <- stats::dist(t(norm_mat))
  dist_mat <- as.matrix(sample_dist)

  plot_df <- as.data.frame(as.table(dist_mat), stringsAsFactors = FALSE)
  names(plot_df) <- c("sample1", "sample2", "distance")

  sample_order <- rownames(dist_mat)[stats::hclust(sample_dist)$order]
  plot_df$sample1 <- factor(plot_df$sample1, levels = sample_order)
  plot_df$sample2 <- factor(plot_df$sample2, levels = rev(sample_order))

  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample1, y = sample2, fill = distance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "steelblue", high = "white") +
    Ol_Reliable() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Euclidean\ndistance")

  plot
}
