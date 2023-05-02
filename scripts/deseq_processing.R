library(fs)
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(gghighlight)
library(ggridges)
library(ggplot2)
library(gghighlight)
library(ggrepel)


bed_read <- function(x){
  read.table(x, header = FALSE, sep="\t", 
             col.names = c("chr", "start", "end", "coverage"), 
             stringsAsFactors=FALSE, quote="")
}

process_expression_data <- function(folder_path, output_file_path) {
  # get a list of all files in the folder that end with .tsv
  expressionfiles <- dir_ls(path = folder_path, regexp = "*.tsv$")

  # read all expression data files and combine them into a single data frame
  expressionData <- expressionfiles |>
    map_df(read_tsv, .id = "filename") |>
    as_tibble() |>
    mutate(
      Filename = as.factor(filename),
      Sample = as.factor(str_extract(filename, "Newman.+(?=\\.bam)")),
    ) |>
    dplyr::rename(Gene = `#Feature`) |>
    select(Sample, Counts, Gene) |>
    mutate(Gene = str_replace(Gene, "cds-|UAKGERNU_", "")) |>
    pivot_wider(names_from = Sample, values_from = Counts)

  # write expression data to file
  write_tsv(expressionData, file = output_file_path)

  return(expressionData)
}

process_dds <- function(metadata_file_path, expression_file_path) {
  # read metadata file
  metadata <- read.table(metadata_file_path, sep = "\t", header = T, row.names = 1)

  # read expression data file
  expressionData <- read.table(expression_file_path, sep = "\t", header = T, row.names = 1, check.names = F)

  # select and order columns in expression data
  expressionDataOrdered <- expressionData |>
    select(
      "Newman-control-R01",
      "Newman-control-R02",
      "Newman-control-R03",
      "Newman-min2-R01",
      "Newman-min2-R02",
      "Newman-min2-R03",
      "Newman-min5-R01",
      "Newman-min5-R02",
      "Newman-min5-R03",
      "Newman-min10-R01",
      "Newman-min10-R02",
      "Newman-min10-R03",
      "Newman-min20-R01",
      "Newman-min20-R02",
      "Newman-min20-R03",
      "Newman-min30-R01",
      "Newman-min30-R02",
      "Newman-min30-R03"
    )

  # check if column names match row names in metadata
  all(colnames(expressionData) == rownames(metadata))

  # convert expression data to a matrix
  expValues <- data.matrix(expressionDataOrdered)

  # round expression data to integers
  expValuesInteger <- round(expValues)

  # create DESeqDataSet object
  ddsObject <- DESeqDataSetFromMatrix(countData = expValuesInteger, colData = metadata, design = ~time)

  # run DESeq analysis
  dds <- DESeq(ddsObject)

  return(dds)
}

get_fold_change_data <- function(dds) {
  time_M2_vs_C <- lfcShrink(dds, coef = "time_M2_vs_C", type = "normal", lfcThreshold = 1) |>
    as_tibble(rownames = "gene") |>
    mutate(sample = "M2")

  time_M5_vs_C <- lfcShrink(dds, coef = "time_M5_vs_C", type = "normal", lfcThreshold = 1) |>
    as_tibble(rownames = "gene") |>
    mutate(sample = "M5")

  time_M10_vs_C <- lfcShrink(dds, coef = "time_M10_vs_C", type = "normal", lfcThreshold = 1) |>
    as_tibble(rownames = "gene") |>
    mutate(sample = "M10")

  time_M20_vs_C <- lfcShrink(dds, coef = "time_M20_vs_C", type = "normal", lfcThreshold = 1) |>
    as_tibble(rownames = "gene") |>
    mutate(sample = "M20")

  time_M30_vs_C <- lfcShrink(dds, coef = "time_M30_vs_C", type = "normal", lfcThreshold = 1) |>
    as_tibble(rownames = "gene") |>
    mutate(sample = "M30")

  foldChangeData <- rbind(time_M2_vs_C, time_M5_vs_C, time_M10_vs_C, time_M20_vs_C, time_M30_vs_C)

  return(foldChangeData)
}

get_fold_change_matrix <- function(dds) {
  time_M2_vs_C <- lfcShrink(dds, coef = "time_M2_vs_C", type = "normal", lfcThreshold = 1)
  time_M5_vs_C <- lfcShrink(dds, coef = "time_M5_vs_C", type = "normal", lfcThreshold = 1)
  time_M10_vs_C <- lfcShrink(dds, coef = "time_M10_vs_C", type = "normal", lfcThreshold = 1)
  time_M20_vs_C <- lfcShrink(dds, coef = "time_M20_vs_C", type = "normal", lfcThreshold = 1)
  time_M30_vs_C <- lfcShrink(dds, coef = "time_M30_vs_C", type = "normal", lfcThreshold = 1)

  M2 <- as_tibble(time_M2_vs_C, rownames = "gene") |>
    mutate(sample = "M2") |>
    pivot_longer(cols = !c(gene, sample), names_to = "parameter", values_to = "value")

  M5 <- as_tibble(time_M5_vs_C, rownames = "gene") |>
    mutate(sample = "M5") |>
    pivot_longer(cols = !c(gene, sample), names_to = "parameter", values_to = "value")

  M10 <- as_tibble(time_M10_vs_C, rownames = "gene") |>
    mutate(sample = "M10") |>
    pivot_longer(cols = !c(gene, sample), names_to = "parameter", values_to = "value")

  M20 <- as_tibble(time_M20_vs_C, rownames = "gene") |>
    mutate(sample = "M20") |>
    pivot_longer(cols = !c(gene, sample), names_to = "parameter", values_to = "value")

  M30 <- as_tibble(time_M30_vs_C, rownames = "gene") |>
    mutate(sample = "M30") |>
    pivot_longer(cols = !c(gene, sample), names_to = "parameter", values_to = "value")

  foldChangeData <- rbind(M2, M5, M10, M20, M30) |>
    filter(parameter == "log2FoldChange") |>
    select(gene, sample, value) |>
    pivot_wider(names_from = sample, values_from = value)

  foldChangeMatrix <- as.matrix(foldChangeData[, -1])

  rownames(foldChangeMatrix) <- foldChangeData$gene

  return(foldChangeMatrix)
}

circular_heatmap <- function(foldChangeMatrix, output_file) {
  circos.clear()
  
  hc <- hclust(dist(foldChangeMatrix), method = "complete")
  
  circos.par(start.degree = 90, gap.degree = 3)
  pdf(output_file, width = 10, height = 10)

  # split <- sample(c("a", "b", "c"), nrow(foldChangeMatrix), replace = TRUE)
  # split <- factor(split, levels = c("a", "b", "c"))

  col_fun1 <- colorRamp2(c(-5, 0, 10), c("#467CBC", "#F8F7CE", "#D22729"))
  lgd <- Legend(title = "Expression", col_fun = col_fun1)

  circos.heatmap(foldChangeMatrix,
    col = col_fun1, split = 3, show.sector.labels = TRUE, rownames.side = "outside",
    dend.side = "inside",
    dend.callback = function(dend, m, si) {
      dendsort(dend)
    },
    track.height = 0.2
  )
  grid.draw(lgd)
  dev.off()
}


volcano_plot <- function(foldChangeData, output_file) {
  
  phanotatefoldChangeDataLog <- foldChangeData %>%
    mutate(log_pValue = -log10(padj)) %>%
    select(gene, log_pValue, log2FoldChange, sample)
  
  ggplot(phanotatefoldChangeDataLog, aes(x = log2FoldChange, y = log_pValue, label = gene)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    geom_vline(xintercept = c(-2,2),linetype = "dotted", color = "red" ) +
    facet_wrap(~factor(sample, levels=c('M2', 'M5', 'M10', 'M20', 'M30')), ncol = 1) +
    annotate("text", label = "P value < 0.05", y = 2, x = 5.5, size = 3, color = "red") +
    gghighlight(log_pValue > -log10(0.05) & abs(log2FoldChange) > 2 , label_key = gene) +
    geom_text_repel(data = subset(phanotatefoldChangeDataLog, log_pValue > -log10(0.05)), segment.size = 0.2, segment.color = "grey50", box.padding = 0.5, min.segment.length = 0) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      strip.background = element_blank(),
      legend.position = "right"
    ) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_y_continuous(breaks = seq(0,6,1)) +
    labs(
      x = "Log2(Fold change)",
      y = "-Log10(Adjusted P value)",
      title = "Volcano plots",
      subtitle = "Differential gene expression of phageK during Newman strain at different sample points",
      caption = ""
    )
  
  ggsave(output_file,  width = 10, height = 12)
}
