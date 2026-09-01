Gliob <- read.csv("Downloads/GBMgenes.csv", stringsAsFactors = FALSE)
Norms <- read.csv("Downloads/normalgenes.csv", stringsAsFactors = FALSE)
str(control_expr)

# We set the genes as row names
gbm_expr <- as.matrix(Gliob[,-1])
control_expr <- as.matrix(Norms[,-1])

rownames(gbm_expr) <- Gliob$Gene
rownames(control_expr) <- Norms$Gene

#The t-test for each Gene
gene_list <- rownames(gbm_expr)

p_values <- sapply(gene_list, function(gene) {
  x <- as.numeric(gbm_expr[gene, ])
  y <- as.numeric(control_expr[gene, ])
  t.test(x, y)$p.value
})

head(p_values)

# We compute the logFC

mean_gbm  <- rowMeans(gbm_expr)
mean_ctrl <- rowMeans(control_expr)

logFC <- log2((mean_gbm + 1) / (mean_ctrl + 1))


head(logFC)


# adjusting p-values with Benjamini–Hochberg

adj_p  <- p.adjust(p_values, method="BH")
head(adj_p)

results <- data.frame(
  Gene     = gene_list,
  logFC    = logFC,
  P.Value  = p_values,
  adj.P.Val= adj_p,
  stringsAsFactors = FALSE
)

head(results)

# sorting to get significance
results <- results[order(results$adj.P.Val), ]
head(results,10)
length(results$Gene)

# Thresholds
padj_threshold <- 0.05
logFC_threshold <- 1

# Filter significant genes
significant_genes <- subset(results, adj.P.Val < padj_threshold &
                              abs(logFC) > logFC_threshold)
head(significant_genes)
length(significant_genes$Gene)

# Checking how many were selected
nrow(significant_genes)

results$Significance <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, 
                               "Significant", "Not Significant")
######################################################
# Volcano plot code

results$negLogP <- -log10(results$adj.P.Val)

ggplot(results, aes(x = logFC, y = negLogP, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(
    x = "Log2FC",
    y = "-Log_10(P-value)",
    color = "Significance") +
  theme_minimal()


library(ggrepel)
ggplot(results, aes(x = logFC, y = negLogP, label = ifelse(adj.P.Val < 0.05 &
                                                             abs(logFC) > 1, Gene, ""))) +
  geom_point(alpha = 0.6) +
  geom_text_repel() +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed") +
  theme_minimal()
##########################################################################
###########################################################################
####################### 350 FOR NEUROIMAGENE ##############################
###########################################################################

next_gene_ <- the_350_genes$Gene

next_gene_ <- as.character(next_gene_)


NextF350 <- neuroimaGene(
  gene_list = next_gene_,
  atlas = "Desikan",
  mtc = 'BH',
  filename="NextF350_results.csv" 
)
head(NextF350)
dim(NextF350)
str(NextF350)

NextF350_nidps <- plot_nidps(NextF350, maxNidps = 20)


length(unique(NextF350$gene_name))
uniq_brainregions <- unique(NextF350$gwas_phenotype)
length(uniq_brainregions)  
head(uniq_brainregions)  

write.csv(uniq_brainregions, "uniq350_region_phenotypes.csv", row.names = FALSE)


################# SIGNIFICANT ASSOCIATION ####################################

sig_assoc <- NextF350[NextF350$atl_BHpval < 0.05, ]

sig_assoc <- sig_assoc[order(abs(sig_assoc$zscore), decreasing = TRUE), ]
dim(sig_assoc)
nrow(sig_assoc)
head(sig_assoc, 15)
tail(sig_assoc$zscore,60)
head(sig_assoc, 10)

top_regions <- unique(sig_assoc$gwas_phenotype)[1:50]
length(unique(top_regions))
top_regions