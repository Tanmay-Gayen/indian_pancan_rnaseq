## Loading packages

library(data.table)
library(ggh4x)
library(ggrepel)
library(ggsci)
library(readxl)
library(stringr)
library(tidyverse)

## Data manipulation for gene-level summaries

# data_expression = fread("data.csv")
# data_grouping = fread("labels.csv", header = T)
# data_genes = read_xlsx("genes.xlsx")
#
# data_combined = left_join(data_grouping, data_expression, by = "V1")
# data_combined = data_combined %>%
#   pivot_longer(cols = 3:20533,
#                names_to = "gene",
#                values_to = "expression")
# colnames(data_combined)[1:2] = c("id", "cancer")
# data_combined = left_join(data_combined, data_genes, by = "gene") %>%
#   mutate(gene = gene_name) %>%
#   select("id", "cancer", "gene", "expression") %>%
#   filter(!(is.na(gene)) &
#            str_detect(gene, "[A-Za-z]"))
#
# save(data_combined, file = "data.rda")

load("data.rda")

data_summary = data_combined %>%
  group_by(cancer, gene) %>%
  summarise(exp_mean = mean(expression, na.rm = T),
            exp_sd = sd(expression, na.rm = T))

data_summary_2 = data_summary %>%
  group_by(gene) %>%
  summarise(exp_mean2 = mean(exp_mean, na.rm = T),
            exp_se = sd(exp_mean, na.rm = T))
data_summary_2$label = ifelse(data_summary_2$exp_se > 5, data_summary_2$gene, "")
data_summary_2$color = ifelse(data_summary_2$label == "", "grey", "black")

## Cross-cancer cross-gene scatterplot

figure_1 = ggplot(data = data_summary_2,
                  mapping = aes(x = exp_mean2, y = exp_se, color = color)) +
  geom_point() +
  theme_light() +
  theme(text = element_text(face = "bold"), legend.position = "none") +
  xlab("Mean Expression Level across Five Cancers") +
  ylab("SE of Means for Five Cancers") +
  geom_text_repel(mapping = aes(label = label), max.overlaps = 100) +
  scale_color_manual(values = c("black", "grey")) +
  ggtitle(label = "Only the genes with an SE of five or higher across the five cancers are labeled.")

if (!(dir.exists("figures"))) {
  dir.create("figures")
}
ggsave(
  filename = "figures/figure_1.pdf",
  plot = figure_1,
  width = 9,
  height = 9
)

## Cross-cancer gene-specific boxplot

important_genes = sort(unique(data_summary_2$label))[-1]
if (!(dir.exists("figures/figure_2"))) {
  dir.create("figures/figure_2")
}

for (specific_gene in important_genes) {
  figure_2 = ggplot(
    data = data_combined %>%
      filter(gene == specific_gene),
    mapping = aes(x = cancer, y = expression)
  ) +
    geom_boxplot() +
    theme_light() +
    theme(text = element_text(face = "bold")) +
    ggtitle(label = paste0("Cross-cancer summary of the expression of ", specific_gene, ".")) +
    xlab("") + ylab("Gene Expression Values for Individual Samples")
  ggsave(
    filename = paste0("figures/figure_2/", specific_gene, ".pdf"),
    plot = figure_2,
    width = 6,
    height = 6
  )
}

## Data manipulation for cross-sample summaries

data_PCA = data_combined %>%
  group_by(id, cancer, gene) %>%
  summarize(expression = max(expression)) %>%
  pivot_wider(values_from = expression, names_from = gene)
data_PCA = data.frame(data_PCA[, 1:2], data_PCA[, data_summary_2$gene[data_summary_2$exp_se > 1]])
PCA = prcomp(as.matrix.data.frame(data_PCA[, -(1:2)]))
data_pcaplot = data.frame(data_PCA[, 1:2], PC1 = PCA$x[, 1], PC2 = PCA$x[, 2])

## Cross-cancer PCA scatterplot

figure_3 = ggplot(data = data_pcaplot, mapping = aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = cancer),
             size = 3,
             alpha = 0.6) +
  theme_light() +
  theme(
    text = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.position = "top"
  ) +
  labs(color = "") +
  xlab(paste0(
    "PC1 (Variance explained: ",
    round(100 * PCA$sdev[1] ^ 2 /
            sum(PCA$sdev ^ 2), 2),
    "%)"
  )) +
  ylab(paste0(
    "PC2 (Variance explained: ",
    round(100 * PCA$sdev[2] ^ 2 /
            sum(PCA$sdev ^ 2), 2),
    "%)"
  )) +
  scale_color_frontiers()
ggsave(
  filename = "figures/figure_3.pdf",
  plot = figure_3,
  width = 9,
  height = 9
)

## Data manipulation for PCA-specific summaries

specific_PC = 1
data_heatmap = data_combined %>%
  filter(gene %in% names(which(abs(PCA$rotation[, specific_PC]) > 0.05)))
data_scaling = data_heatmap %>%
  group_by(gene) %>%
  summarize(
    cutoff_q1 = quantile(expression, 0.25),
    cutoff_q2 = median(expression),
    cutoff_q3 = quantile(expression, 0.75)
  )
data_heatmap = left_join(data_heatmap, data_scaling, by = "gene")
data_heatmap = data_heatmap %>%
  mutate(grading = ifelse(
    expression < cutoff_q1,
    "[Min, Q1)",
    ifelse(
      expression < cutoff_q2 &
        expression >= cutoff_q1,
      "[Q1, Median)",
      ifelse(
        expression < cutoff_q3 &
          expression >= cutoff_q2,
        "[Median, Q3)",
        "[Q3, Max]"
      )
    )
  ))
data_heatmap$grading = factor(data_heatmap$grading, sort(unique(data_heatmap$grading))[c(2:3, 1, 4)])
data_scaling = data_heatmap %>%
  dplyr::filter(cancer == "KIRC") %>%
  group_by(gene) %>%
  summarize(highs = mean(grading == "[Q3, Max]")) %>%
  arrange(-highs)
data_heatmap$gene = factor(data_heatmap$gene, levels = data_scaling$gene)

## Cross-cancer PCA-specific heatmap

figure_4 = ggplot(data = data_heatmap, mapping = aes(x = gene, y = id)) +
  geom_tile(mapping = aes(fill = grading), alpha = 0.72) +
  facet_nested(
    cancer ~ .,
    strip = strip_nested(clip = "off", background_y = element_rect(color = "black")),
    switch = "both",
    scales = "free",
    space = "free"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90),
    text = element_text(face = "bold"),
    strip.placement = "outside",
    legend.position = "bottom",
    strip.text.y.left = element_text(
      size = 6,
      angle = 0,
      margin = margin()
    ),
    strip.text.x.bottom = element_text(
      size = 6,
      angle = 90,
      margin = margin()
    ),
    strip.switch.pad.grid = unit(0, "cm")
  ) +
  xlab("") + ylab("") +
  labs(fill = "Gene-specific Expression") +
  scale_fill_manual(values = c("blue", "seagreen", "darkgreen", "red")) +
  ggtitle(label = paste0("Top Genes Driving Principal Component ", specific_PC, "."))
if (!(dir.exists("figures/figure_4"))) {
  dir.create("figures/figure_4")
}
ggsave(
  filename = paste0("figures/figure_4/PC", specific_PC, ".pdf"),
  plot = figure_4,
  width = 9,
  height = 9
)

## Data manipulation for PCA-specific summaries

specific_PC = 2
data_heatmap = data_combined %>%
  filter(gene %in% names(which(abs(PCA$rotation[, specific_PC]) > 0.05)))
data_scaling = data_heatmap %>%
  group_by(gene) %>%
  summarize(
    cutoff_q1 = quantile(expression, 0.25),
    cutoff_q2 = median(expression),
    cutoff_q3 = quantile(expression, 0.75)
  )
data_heatmap = left_join(data_heatmap, data_scaling, by = "gene")
data_heatmap = data_heatmap %>%
  mutate(grading = ifelse(
    expression < cutoff_q1,
    "[Min, Q1)",
    ifelse(
      expression < cutoff_q2 &
        expression >= cutoff_q1,
      "[Q1, Median)",
      ifelse(
        expression < cutoff_q3 &
          expression >= cutoff_q2,
        "[Median, Q3)",
        "[Q3, Max]"
      )
    )
  ))
data_heatmap$grading = factor(data_heatmap$grading, sort(unique(data_heatmap$grading))[c(2:3, 1, 4)])
data_scaling = data_heatmap %>%
  dplyr::filter(cancer == "PRAD") %>%
  group_by(gene) %>%
  summarize(highs = mean(grading == "[Q3, Max]")) %>%
  arrange(-highs)
data_heatmap$gene = factor(data_heatmap$gene, levels = data_scaling$gene)

## Cross-cancer PCA-specific heatmap

figure_4 = ggplot(data = data_heatmap, mapping = aes(x = gene, y = id)) +
  geom_tile(mapping = aes(fill = grading), alpha = 0.72) +
  facet_nested(
    cancer ~ .,
    strip = strip_nested(clip = "off", background_y = element_rect(color = "black")),
    switch = "both",
    scales = "free",
    space = "free"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90),
    text = element_text(face = "bold"),
    strip.placement = "outside",
    legend.position = "bottom",
    strip.text.y.left = element_text(
      size = 6,
      angle = 0,
      margin = margin()
    ),
    strip.text.x.bottom = element_text(
      size = 6,
      angle = 90,
      margin = margin()
    ),
    strip.switch.pad.grid = unit(0, "cm")
  ) +
  xlab("") + ylab("") +
  labs(fill = "Gene-specific Expression") +
  scale_fill_manual(values = c("blue", "seagreen", "darkgreen", "red")) +
  ggtitle(label = paste0("Top Genes Driving Principal Component ", specific_PC, "."))
if (!(dir.exists("figures/figure_4"))) {
  dir.create("figures/figure_4")
}
ggsave(
  filename = paste0("figures/figure_4/PC", specific_PC, ".pdf"),
  plot = figure_4,
  width = 9,
  height = 9
)