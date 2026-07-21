#' Plot top enrichment results
#'
#' Plots a horizontal bar chart of the top enriched terms from
#' [WebGestalt_Pipe2()] output (ORA or GSEA), ranked by FDR and coloured by
#' direction (enriched vs. depleted).
#'
#' @param enrichment A results data.frame from [WebGestalt_Pipe2()]. Must
#' contain an "ID" (term label) column, an "FDR" column, and either an
#' "enrichmentRatio" column (ORA) or an "NES" column (GSEA).
#' @param n Number of top terms to plot (ranked by FDR, after the `alpha`
#' filter).
#' @param valueCol Column to use for bar length/direction. Defaults to
#' auto-detecting "enrichmentRatio" (ORA) or "NES" (GSEA).
#' @param alpha FDR cut-off; terms above this are excluded before ranking.
#' @return A ggplot2 bar plot object.
#' @importFrom rlang .data
#' @export
plotEnrichment <- function(enrichment, n = 20, valueCol = NULL, alpha = 1) {

  if (is.null(valueCol)) {
    valueCol <- if ("enrichmentRatio" %in% colnames(enrichment)) {
      "enrichmentRatio"
    } else if ("NES" %in% colnames(enrichment)) {
      "NES"
    } else {
      stop("Could not auto-detect a value column. Pass `valueCol` ",
           "explicitly (e.g. \"enrichmentRatio\" or \"NES\").")
    }
  }

  if (!all(c("ID", "FDR", valueCol) %in% colnames(enrichment))) {
    stop("`enrichment` must contain columns \"ID\", \"FDR\", and \"",
         valueCol, "\".")
  }

  plot_df <- enrichment[enrichment$FDR <= alpha, , drop = FALSE]
  plot_df <- plot_df[order(plot_df$FDR), ]
  plot_df <- utils::head(plot_df, n)

  if (nrow(plot_df) == 0) {
    stop("No enrichment terms remain after applying `alpha = ", alpha, "`.")
  }

  plot_df$direction <- ifelse(plot_df[[valueCol]] > 0, "Enriched", "Depleted")
  plot_df$ID <- factor(plot_df$ID, levels = rev(plot_df$ID))

  plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[["ID"]], y = .data[[valueCol]], fill = .data[["direction"]])
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    Ol_Reliable() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
    ggplot2::labs(x = NULL, y = valueCol)

  plot
}
