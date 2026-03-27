## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 7,
  fig.height = 5
)

## ---- message=FALSE-----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

## -----------------------------------------------------------------------------
data("gwas_results", package = "miamiplot")

## -----------------------------------------------------------------------------
gwas_results$consequences <- 
  rep(sample(c("synonymous_variant", "missense_variant", 
               "downstream_gene_variant", "nc_transcript_variant", 
               "upstream_gene_variant", "non_coding_expn_variant", 
               "intron_variant", "NMD_transcript_variant", 
               "3_prime_UTR_variant", NA_character_), nrow(gwas_results)/2, 
             replace = TRUE, prob = c(0.2, 0.18, 0.15, 0.11, 0.1, 0.08, 0.06, 
                                      0.05, 0.03, 0.04)), 2)

## -----------------------------------------------------------------------------
str(gwas_results)

## -----------------------------------------------------------------------------
plot_df <- gwas_results %>% 
  group_by(chr) %>% 
  # Compute chromosome size
  summarise(chrlength = max(pos)) %>%  
  # Calculate cumulative position of each chromosome
  mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
  select(-chrlength) %>%
  # Temporarily add the cumulative length of each chromosome to the initial 
  # dataset 
  left_join(gwas_results, ., by=c("chr"="chr")) %>%
  # Sort by chr then position 
  arrange(chr, pos) %>%
  # Add the position to the cumulative chromosome length to get the position of 
  # this probe relative to all other probes
  mutate(rel_pos = pos + cumulativechrlength) %>%
  # Calculate the logged p-value too
  mutate(logged_p = -log10(pval)) %>%
  select(-cumulativechrlength)

## -----------------------------------------------------------------------------
axis_df <- plot_df %>% 
  group_by(chr) %>% 
  summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)

## -----------------------------------------------------------------------------
maxp <- ceiling(max(plot_df$logged_p, na.rm = TRUE))

## -----------------------------------------------------------------------------
upper_labels<- plot_df %>%
  filter(beta > 0) %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(logged_p)) %>%
  slice(1:5) %>%
  select(rsid)

lower_labels <- plot_df %>%
  filter(beta <= 0) %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(logged_p)) %>%
  slice(1:5) %>%
  select(rsid)

label_list <- rbind(upper_labels, lower_labels)

plot_df <- plot_df %>%
  mutate(label = case_when(rsid %in% label_list$rsid ~ 
                             paste0(rsid, "\n", consequences), TRUE ~ ""))

## -----------------------------------------------------------------------------
upper_plot <- ggplot() + 
  geom_point(data = plot_df[which(plot_df$beta > 0 ),], 
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
             size = 0.25) +
  scale_color_manual(values = rep(c("black", "grey"), nrow(axis_df))) +
  scale_x_continuous(labels = axis_df$chr, 
                     breaks = axis_df$chr_center, 
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, maxp), 
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(1e-5), color = "red", linetype = "dashed", 
             size = 0.3) +
  geom_hline(yintercept = -log10(5e-8), color = "blue", linetype = "solid", 
             size = 0.3) +
  labs(x = "", y = bquote(atop('Positive beta values', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.margin = margin(b = 0, l = 10)) + 
  geom_label_repel(data = plot_df[which(plot_df$beta > 0),],
                   aes(x = rel_pos, y = logged_p, label = label),
                   size = 2, segment.size = 0.2, point.padding = 0.3,
                   ylim = c(maxp / 3, NA), box.padding = 0.5)

upper_plot

## -----------------------------------------------------------------------------
lower_plot <- ggplot() + 
  geom_point(data = plot_df[which(plot_df$beta <= 0),], 
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
             size = 0.25) +
  scale_color_manual(values = rep(c("black", "grey"), nrow(axis_df))) +
  scale_x_continuous(breaks = axis_df$chr_center, 
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits = c(maxp, 0), 
                     expand = expansion(mult = c(0, 0.02))) + 
  geom_hline(yintercept = -log10(1e-5), color = "red", linetype = "dashed", 
             size = 0.3) +
  geom_hline(yintercept = -log10(5e-8), color = "blue", linetype = "solid", 
             size = 0.3) +
  labs(x = "", y = bquote(atop('Negative beta values', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 0, l = 10)) + 
  geom_label_repel(data = 
                     plot_df[which(plot_df$beta <= 0),],
                   aes(x = rel_pos, y = logged_p, label = label),
                   size = 2, segment.size = 0.2, point.padding = 0.3,
                   ylim = c(NA, -(maxp / 3)), box.padding = 0.5)

lower_plot

## -----------------------------------------------------------------------------
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
p

