# =========================
# LIBRARIEs
# =========================
library(data.table)
library(tidyverse)
library(limpa)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(UniProt.ws)
library(dplyr)
library(clusterProfiler)

# =========================
# LOAD DATA (PROTEIN LEVEL)
# =========================
# load data
pr <- fread("D:/PLUS/bEV_proteomics_hpyl_BEV/raw_data/20260305_report.pr_matrix (1).tsv")

# select intensity columns
expr_cols <- grep("Hpyl|BEV", names(pr), value = TRUE)

# create peptide matrix
y.peptide <- as.matrix(pr[, ..expr_cols])

# unique peptide IDs (important!)
rownames(y.peptide) <- paste0(pr$Precursor.Id)

# peptide → protein mapping
protein.id <- pr$Protein.Group







####checkifnormalised


nor <- as.data.frame(y.peptide) |>
  pivot_longer(cols = everything(), names_to = "sample", values_to = "intensity")

a <- ggplot(nor, aes(x = intensity, color = sample)) +
  geom_density(na.rm = TRUE) +
  scale_x_log10() +
  theme_minimal()

ggsave(
  "D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/normaliseddist.png",
  a, width = 10, height = 6, dpi = 300
)

# =========================
# METADATA
# =========================
metadata <- data.frame(sample = colnames(y.peptide))

metadata$condition <- ifelse(grepl("Hpyl", metadata$sample), "Hpyl", "BEV")
metadata$batch <- sub(".*_(B[0-9]+)_.*", "\\1", metadata$sample)
metadata$replicate <- as.numeric(sub(".*_Repl([0-9]+).*", "\\1", metadata$sample))


metadata <- metadata |>                      
  mutate(condition = fct_relevel(condition, "Hpyl"))


# =========================
 ##LOG2 TRANSFORM
# =========================
expr_log <- log2(y.peptide)



cor_mat <- cor(expr_log, use = "pairwise.complete.obs", method = "pearson")

png("D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/heatmapcorMT.png",
    width = 1600, height = 1600, res = 300)

# Plot the heatmap
pheatmap(
  cor_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Sample Correlation Heatmap"
)

dev.off()


#####plot average missing values

png("D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/avgmissingvalues.png",
    width = 1600, height = 1600, res = 300)

plotAveVsMis(expr_log)

dev.off()



###Keeping only rows rows where atleast 2 values are present

# logical indices for each condition
bev_cols  <- metadata$sample[metadata$condition == "BEV"]
hpyl_cols <- metadata$sample[metadata$condition == "Hpyl"]

# ensure column order matches expr_log
bev_idx  <- colnames(expr_log) %in% bev_cols
hpyl_idx <- colnames(expr_log) %in% hpyl_cols

# apply filtering rule
keep <- rowSums(!is.na(expr_log[, bev_idx])) >= 2 &
  rowSums(!is.na(expr_log[, hpyl_idx])) >= 2

# subset data
expr_log <- expr_log[keep, ]
protein.id <- protein.id[keep]

####Detection probablity curve--------
dpcest <- dpcCN(expr_log)
dpcest$dpc

png("D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/dpccurve.png",
    width = 1600, height = 1600, res = 300)

plotDPC(dpcest)

dev.off()


y <- dpcQuant(
  expr_log,
  protein = protein.id,
  dpc = dpcest
)
dim(y)



png("D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/yflit.png",
    width = 2500, height = 2000, res = 300)

plotMDSUsingSEs(y)

dev.off()




design <- model.matrix(~condition, data = metadata)
fit <- dpcDE(y, design, plot=TRUE)

fit <- eBayes(fit)
summary(dt <- decideTests(fit[,2]))




results <- topTable(
  fit,
  coef = "conditionBEV",
  number = Inf,
  sort.by = "P"
) |>
  rownames_to_column("Protein_ID") |>
  arrange(adj.P.Val)

up <- UniProt.ws(taxId = 210)  # 210 = Helicobacter pylori
ids <- results$Protein_ID
mapped <- AnnotationDbi::select(
  up,
  keys = ids,
  columns = c("accession", "gene_primary", "protein_name"),
  keytype = "UniProtKB"
)

results$Gene <- sapply(results$Protein, function(x) {
  ids <- unlist(strsplit(x, ";"))
  
  genes <- mapped$Gene.Names..primary.[mapped$From %in% ids]
  
  genes <- unique(genes[!is.na(genes)])
  
  if(length(genes) == 0) NA else paste(genes, collapse = ";")
})



volcano_df <- results |>
  mutate(
    negLogP = -log10(adj.P.Val),
    sig = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

bev_labels <- volcano_df |>
  filter(adj.P.Val < 0.05, logFC > 1) |>
  arrange(adj.P.Val, desc(logFC)) |>
  slice_head(n = 40)


p1 <- ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = bev_labels,
    aes(label = Protein),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot (BEV vs Hpyl)",
    x = "Log2 Fold Change",
    y = "-Log10 adj.P value",
    color = "Significance"
  )

ggsave(
  "D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/volcano_proteins.png",
  p1, width = 10, height = 6, dpi = 300
)


b <- ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # ---- BEV-only labels ----
geom_text_repel(
  data = bev_labels,
  aes(label = Gene),
  size = 3,
  max.overlaps = Inf
) +
  theme_minimal() +
  labs(
    title = "Volcano Plot (BEV vs Hpyl)",
    x = "Log2 Fold Change",
    y = "-Log10 adj.P value",
    color = "Significance"
  )

ggsave("D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/volcanoDEGlimpagene.png",
       plot = b, width = 10, height = 6, dpi = 300)





######Enrichmentanalysis

write.table(
  bev_labels$Gene,
  "bev_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
writeLines(bev_labels$Gene, "D:/PLUS/bEV_proteomics_hpyl_BEV/results/limpa/bev_genes.txt")
