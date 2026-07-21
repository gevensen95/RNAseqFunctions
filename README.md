# RNAseqBasic

Convenience functions for RNA-seq differential expression workflows: run DESeq2/edgeR from raw counts, visualize results, convert gene IDs, and run WebGestalt enrichment — all in a few consistent function calls.

## Installation

```r
# install.packages("devtools")
devtools::install_local("gevensen95/RNAseqBasic")
```

RNAseqBasic itself only requires CRAN/Bioconductor packages available through normal installation (`biomaRt`, `ggplot2`, `ggrepel`, `rlang`, `WebGestaltR`). The actual differential expression engines are optional and installed separately depending on which one you use:

```r
# install.packages("BiocManager")
BiocManager::install("DESeq2")   # for method = "DESeq2" (default)
BiocManager::install("edgeR")    # for method = "edgeR"
BiocManager::install("limma")    # only needed for edgeR + a custom `contrast`
BiocManager::install("SummarizedExperiment")  # used internally alongside DESeq2
```

You only need to install the engine(s) you actually plan to use — `makeComp()` and friends will tell you which package is missing if you try to use one you don't have.

## Tutorial

This walks through a full workflow on simulated count data: QC, differential expression, visualization, gene ID conversion, and enrichment.

### 1. Set up example data

`makeComp()` and the plotting functions expect a raw counts matrix (genes × samples) and a metadata data.frame whose row names match the counts' column names exactly.

```r
library(RNAseqBasic)

set.seed(42)
n_genes <- 2000
n_samples <- 10

counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 200, size = 1),
  nrow = n_genes,
  dimnames = list(paste0("gene", seq_len(n_genes)), paste0("sample", seq_len(n_samples)))
)

metadata <- data.frame(
  condition = factor(rep(c("control", "treated"), each = 5),
                      levels = c("control", "treated")),
  row.names = paste0("sample", seq_len(n_samples))
)
```

### 2. QC before you test anything

`PCAplot()` and `sampleDistHeatmap()` both normalize the counts for you (DESeq2's `vst()` or edgeR's log-CPM) and are worth a look before running any comparison, to catch batch effects or outlier samples early.

```r
PCAplot(counts, metadata, color_by = "condition")
sampleDistHeatmap(counts, metadata)
```

### 3. Run the differential expression comparison

`makeComp()` wraps either DESeq2 (default) or edgeR behind one interface. If `metadata` has a `condition` column, the `~condition` design is used automatically — pass `formula` explicitly for anything more complex.

```r
res <- makeComp(counts, metadata)
head(res)
```

For a specific contrast, or to use edgeR instead:

```r
# DESeq2 with an explicit contrast
res <- makeComp(counts, metadata, contrast = c("condition", "treated", "control"))

# edgeR instead of DESeq2
res_edger <- makeComp(counts, metadata, method = "edgeR")
```

Every `makeComp()` result is shaped the same way regardless of engine — `log2FoldChange` in column 2, gene IDs in a "Gene" column in column 7 — so it plugs directly into every function below.

### 4. Get a quick summary

```r
summarizeDEGs(res, alpha = 0.05, FCsig = 1)
#>   comparison up down not_significant total
#> 1 comparison 84   79            1837  2000
```

Running several comparisons? Pass a named list and get one row per comparison:

```r
summarizeDEGs(list(treated_vs_control = res1, timepoint2_vs_timepoint1 = res2))
```

### 5. Visualize the results

```r
eRupt2(res, alpha = 0.05, FCsig = 1)                              # volcano plot
MAplot(res, alpha = 0.05, FCsig = 1)                               # mean expression vs. fold change
DEGheatmap(counts, metadata, res, n = 30, annotate_by = "condition")  # top-30 DEG heatmap
```

All plots share the same look via the package's built-in `Ol_Reliable()` ggplot2 theme, so figures stay consistent across a project without any extra styling code.

### 6. Convert gene IDs

`convert_gene_names()` and `convert_ensembl_to_symbol()` query Ensembl BioMart (requires internet access) to map between IDs, symbols, and other attributes for mouse or human.

```r
convert_ensembl_to_symbol(c("ENSMUSG00000064341", "ENSMUSG00000051951"), species = "mouse")
```

`convert_gene_names()` is the more general form — same idea, but lets you set the filter type (e.g. look up by `"mgi_symbol"` instead of Ensembl ID) and request extra BioMart attributes:

```r
convert_gene_names(c("Actb", "Gapdh"), species = "mouse", gene_filters = "mgi_symbol",
                    extra_attributes = "chromosome_name")
```

### 7. Run enrichment analysis

`WebGestalt_Pipe2()` looks genes up **by gene symbol**, so if your counts (and therefore your `makeComp()` results) are indexed by Ensembl ID, convert the "Gene" column to symbols first:

```r
symbols <- convert_ensembl_to_symbol(res$Gene, species = "mouse")
res$Gene <- symbols$mgi_symbol[match(res$Gene, symbols$ensembl_gene_id)]
res <- res[!is.na(res$Gene) & res$Gene != "", ]

enrich <- WebGestalt_Pipe2(mode = "ORA", results = res, alpha = 0.05, FC = 1)
plotEnrichment(enrich)
```

Switch to `mode = "GSEA"` to rank on log2FoldChange directly instead of thresholding into up/down gene sets:

```r
gsea <- WebGestalt_Pipe2(mode = "GSEA", results = res, DB = "pathway_KEGG")
plotEnrichment(gsea)   # auto-detects the "NES" column for GSEA results
```

Both `WebGestalt_Pipe2()` and `convert_gene_names()`/`convert_ensembl_to_symbol()` make live requests to external services (WebGestalt and Ensembl BioMart, respectively) — they need an internet connection and can be slow on large gene lists.

## Function reference

**QC**

| Function | What it does |
|---|---|
| `PCAplot()` | Sample PCA from normalized counts, colour/shape by metadata |
| `sampleDistHeatmap()` | Hierarchically clustered sample-to-sample distance heatmap |

**Differential expression**

| Function | What it does |
|---|---|
| `makeComp()` | Run DESeq2 or edgeR from raw counts + metadata + a design formula |
| `summarizeDEGs()` | Tally up/down/not-significant genes across one or more comparisons |

**Visualization**

| Function | What it does |
|---|---|
| `eRupt2()` | Volcano plot |
| `MAplot()` | Mean expression vs. log2 fold change |
| `DEGheatmap()` | Z-scored expression heatmap of the top DEGs |
| `Ol_Reliable()` | The shared ggplot2 theme used across all plots above |

**Gene ID conversion**

| Function | What it does |
|---|---|
| `convert_gene_names()` | General Ensembl BioMart ID/attribute lookup |
| `convert_ensembl_to_symbol()` | Shortcut: Ensembl ID → gene symbol |

**Enrichment**

| Function | What it does |
|---|---|
| `WebGestalt_Pipe2()` | ORA or GSEA enrichment via WebGestaltR |
| `plotEnrichment()` | Bar plot of the top enriched terms |

## License

MIT — see [LICENSE.md](LICENSE.md).
