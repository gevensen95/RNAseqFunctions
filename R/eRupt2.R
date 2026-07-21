#' Plots a Volcano Plot
#'
#' This function plots a volcano plot for RNAseq differential expression
#' results, colouring points by significance and labelling the top 10 hits.
#'
#' @param results Dataframe of results. Must contain columns named "padj"
#' and "log2FoldChange", plus a gene identifier column (see `geneID`).
#' @param alpha alpha value for significance
#' @param FCsig Fold change value for significance
#' @param geneID Column for labeling genes
#' @return A ggplot2 volcano plot object
#' @export

eRupt2 <- function(results, alpha = 0.05, FCsig = 1, geneID = 'Gene'){

  #colour significance by p-adj:
  results$significance <- "NS"
  results$significance[results$padj < alpha & results$log2FoldChange > FCsig] <- "Upregulated"
  results$significance[results$padj < alpha & results$log2FoldChange < (FCsig*-1)] <- "Downregulated"
  results <- stats::na.omit(results)
  results$significance <- factor(results$significance, levels=c("Downregulated", "NS", "Upregulated"))

  # changes FDR=0 to very small values based on the rest of the data, so log transformed FDR=0 values are still plotted
  min_val <- min(stats::na.omit(results$padj[results$padj != 0]))
  results$padj[results$padj == 0] <- stats::runif(length(stats::na.omit(results$padj[results$padj == 0])), (min_val*1e-3), (min_val*1e-1))
  results <- results[order(-log10(results$padj), decreasing = TRUE), ]
  results$label <- FALSE
  results$label[seq_len(min(10, nrow(results)))] <- TRUE

  plot <- ggplot2::ggplot(data = results, ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = significance))+
    ggplot2::geom_point(alpha=0.35)+
    Ol_Reliable()+
    ggplot2::scale_y_continuous( limits=c(0, (-log10(min_val)+(-log10(min_val)*0.1))) )+
    ggplot2::theme(legend.title = ggplot2::element_blank())+
    ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(label, results[[geneID]], "")), size=2, max.overlaps = Inf)+
    ggplot2::scale_color_manual(values=c("blue","black","red"), drop=F)

  return(plot)

}
