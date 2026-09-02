###### NORMAL HUMAN BRAIN GENE DATA ############### GSE50161 14 samples ##########
normal_data <- read.delim("~/Downloads/GSE50161_series_matrix.txt", comment.char = "!", header = TRUE, sep = "\t")

head(normal_data[,1:15], 15)

probe_ids <- normal_data[[1]] 

nomgene_symb <- mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

head(nomgene_symb, 10)

normal_data$Gene <- nomgene_symb

head(normal_data[, c("Gene", "GSM1214845", "GSM1214843", "GSM1214836", "GSM1214838")], 10)

#Creating my Normal Genes Dataframe

normal_samples <- c("GSM1214936", "GSM1214937", "GSM1214938", "GSM1214939", "GSM1214940","GSM1214941",
                    "GSM1214942", "GSM1214943", "GSM1214944", "GSM1214945", "GSM1214946","GSM1214947",
                    "GSM1214948")


My_norm_genes <- normal_data[, c("Gene", normal_samples)]

head(My_norm_genes)

# Convert all columns except 'Gene' to numeric
My_norm_genes[, -1] <- lapply(My_norm_genes[, -1], as.numeric)

My_norm_genes

str(My_norm_genes)
My_norm_genes$Gene
write.csv(My_norm_genes, "Normal_Genes.csv")

################MERGING########################

#common_genes <- intersect(My_norm_genes$Gene, gbm_45samples$Gene)

#My_norm_genes_aligned <- My_norm_genes[match(common_genes, My_norm_genes$Gene), ]
#gbm_genes_aligned <- gbm_45samples[match(common_genes, gbm_45samples$Gene), ]
#head(My_norm_genes_aligned)
#head(gbm_genes_aligned)

