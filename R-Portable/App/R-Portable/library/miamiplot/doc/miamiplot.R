## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 7,
  fig.height = 5
)

## ---- message=FALSE-----------------------------------------------------------
library(miamiplot)
library(dplyr)
library(ggplot2)
requireNamespace("RColorBrewer", quietly = TRUE)

## -----------------------------------------------------------------------------
str(gwas_results)

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results[which(gwas_results$study == "A"),], 
        split_by = "beta", split_at = 0, p = "pval",
        upper_ylab = "Positive beta values",
        lower_ylab = "Negative beta values")

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B")

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A",lower_ylab = "Study B", suggestive_line = NULL, 
        genome_line = 5e-15)

## -----------------------------------------------------------------------------
mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, 
                                                                 "Paired"))(22)

ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", chr_colors = mycolors, 
        genome_line_color = "black", suggestive_line_color = "#A9A9A9")

## -----------------------------------------------------------------------------
my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[3:4]

ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", chr_colors = NULL,
        upper_chr_colors = my_upper_colors, lower_chr_colors = my_lower_colors,
        genome_line_color = "black", suggestive_line_color = "#A9A9A9")

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", 
        hits_label_col = "rsid")

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", 
        hits_label_col = c("rsid", "beta"), top_n_hits = 10)

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", 
        hits_label = c("rs19142", "rs27017", "rs19103", "rs26991", "rs29240"),
        hits_label_col = "rsid")

## -----------------------------------------------------------------------------
plot_data <- prep_miami_data(data = gwas_results, split_by = "study", 
                             split_at = "A", p = "pval")

## -----------------------------------------------------------------------------
# Study A
studyA_labels <- plot_data$upper %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0(rsid, "\n", beta)) %>%
  select(rel_pos, logged_p, label) %>%
  arrange(desc(logged_p)) %>%
  slice(1:5)

# Study B
studyB_labels <- plot_data$lower %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0(rsid, "\n", beta)) %>%
  select(rel_pos, logged_p, label) %>%
  arrange(desc(logged_p)) %>%
  slice(1:5)

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", 
        upper_labels_df = studyA_labels, lower_labels_df = studyB_labels)

## -----------------------------------------------------------------------------
# Get the position of the two peaks +- 100 bp.
studyA_highlight_pos <- plot_data$upper %>%
  filter(pval < 5e-8) %>%
  group_by(chr) %>%
  filter(pval == min(pval)) %>%
  summarise(start = rel_pos - 100, end = rel_pos + 100) %>%
  select(-chr) %>%
  apply(., 1, function(x) x["start"]:x["end"]) %>%
  as.vector()

# Find which rsids match these SNPs
studyA_highlight_rsid <- plot_data$upper %>%
  mutate(in_peak = case_when(rel_pos %in% studyA_highlight_pos ~ "Yes", 
                             TRUE ~ "No")) %>%
  filter(in_peak == "Yes") %>%
  select(rsid)

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", 
        upper_labels_df = studyA_labels, lower_labels_df = studyB_labels, 
        upper_highlight = studyA_highlight_rsid$rsid, 
        upper_highlight_col = "rsid")

## -----------------------------------------------------------------------------
studyA_highlight_rsid <- plot_data$upper %>%
  mutate(in_peak = case_when(rel_pos %in% studyA_highlight_pos ~ "Yes", 
                             TRUE ~ "No")) %>%
  filter(in_peak == "Yes") %>%
  arrange(logged_p) %>%
  mutate(color = rep(c("magenta", "green"), length.out = n())) %>%
  select(rsid, logged_p, color)

## -----------------------------------------------------------------------------
ggmiami(data = gwas_results, split_by = "study", split_at = "A", p = "pval", 
        upper_ylab = "Study A", lower_ylab = "Study B", 
        upper_labels_df = studyA_labels, lower_labels_df = studyB_labels, 
        upper_highlight = studyA_highlight_rsid$rsid, 
        upper_highlight_col = "rsid", 
        upper_highlight_color = studyA_highlight_rsid$color)

## -----------------------------------------------------------------------------
p <- ggmiami(data = gwas_results, split_by = "study", split_at = "A", 
             p = "pval", upper_ylab = "Study A", lower_ylab = "Study B", 
             upper_labels_df = studyA_labels, lower_labels_df = studyB_labels, 
             upper_highlight = studyA_highlight_rsid$rsid, 
             upper_highlight_col = "rsid", upper_highlight_color = "magenta")

# ggsave(p, filename = "ExampleMiamiPlot.png", device = "png", width = 8, 
       # height = 6, units = "in")

