####################################################################################
########## GSE36245 46 samples #################################################

expr_data <- read.delim("~/Downloads/GSE36245_series_matrix.txt", comment.char = "!", header = TRUE, sep = "\t")

head(expr_data[,1:5], 15)

library(BiocManager)

BiocManager::install("hgu133plus2.db")

library(hgu133plus2.db)  #used to map probe IDs to gene symbols

gene_symbols <- mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

head(gene_symbols)

expr_data$Gene <- gene_symbols

head(expr_data[, c("Gene", "GSM884997", "GSM884998", "GSM884999", "GSM885000")], 10)

head(expr_data)

#keytypes(hgu133plus2.db)

