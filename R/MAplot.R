#' Plots an MA plot
#'
#' Plots mean expression against log2 fold change for a DEG results
#' data.frame (e.g. from [makeComp()]), colouring points by significance.
#' Significant points are drawn on top of the non-significant background.
#'
#' @param results Dataframe of DEG results. Must contain columns named
#' "baseMean", "log2FoldChange", and "padj" (as returned by [makeComp()]).
#' @param alpha padj cut-off for significance
#' @param FCsig absolute log2FC cut-off for significance
#' @return A ggplot2 MA plot object.
#' @export
MAplot <- function(results, alpha = 0.05, FCsig = 0) {

  results$significance <- "NS"
  results$significance[results$padj < alpha & results$log2FoldChange > FCsig] <- "Upregulated"
  results$significance[results$padj < alpha & results$log2FoldChange < (FCsig * -1)] <- "Downregulated"
  results <- stats::na.omit(results)
  results$significance <- factor(results$significance, levels = c("Downregulated", "NS", "Upregulated"))

  # draw significant points last so they sit on top of the NS cloud
  results <- results[order(results$significance == "NS", decreasing = TRUE), ]

  plot <- ggplot2::ggplot(results, ggplot2::aes(x = baseMean, y = log2FoldChange, color = significance)) +
    ggplot2::geom_point(alpha = 0.35) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    Ol_Reliable() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("blue", "black", "red"), drop = FALSE)

  plot
}
