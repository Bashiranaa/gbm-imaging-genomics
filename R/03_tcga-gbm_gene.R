############################################################################
################## GENE DATA PREPROCESSING ######################################
############################################################################
###### TCGA-GBM DATA GSE83130 #############


# Read and tidy each file
expression_list <- map(files, function(file) {
  df <- read_tsv(file, skip = 1, col_names = c("Gene", "Expression"))
  sample_id <- str_remove(file, "_.*")  # extract GSM ID
  df <- df %>% mutate(Sample = sample_id)
})

str(expression_list, max.level = 1)
expression_list[[1]]

# Combine all tibbles into one big long tibble
combined_expr <- bind_rows(expression_list)

# Remove any header rows that snuck in
combined_expr <- combined_expr %>%
  filter(Gene != "Composite Element REF")

# Convert expression to numeric
combined_expr <- combined_expr %>%
  mutate(Expression = as.numeric(Expression))

# Pivot to wide format
expression_matrix <- combined_expr %>%
  pivot_wider(names_from = Sample, values_from = Expression)

dim(expression_matrix)      # Rows = genes, Columns = samples (+1 for gene name)

expression_matrix[1:10, 1:15] # Peek at the first few genes and samples

write_csv(expression_matrix, "expression_matrix.csv")

boxplot(expression_matrix[,-1], main = "Raw Expression Distributions", las = 2)

head(expression_matrix)
dim(expression_matrix)

selected_indices <- c(15, 16, 19, 20, 52, 53, 54, 55, 95, 103, 138, 141, 170, 209,
                      235, 259, 260, 262, 286, 297, 298, 333, 335, 336, 339, 361,
                      377, 388, 389, 340, 408, 409, 410, 440, 444, 455, 463, 466,
                      467, 497, 501, 502, 518, 568, 582)

selected_columns <- c(1, selected_indices)

# Create the new expression matrix
gbm_45samples <- expression_matrix[, selected_columns]

length(gbm_45samples)

str(gbm_45samples)

write.csv(gbm_45samples, "GBM_Genes.csv")
#############################################################################################
###################TESTING FOR SIGNIFICANTLY EXPRESSED GENES ################################
#############################################################################################


# Extract sample names (skip the first column which is "Gene")
sample_ids <- colnames(expression_matrix)[-1]

# Create data frame with a blank column for condition labels
sample_info <- tibble(
  SampleID = sample_ids,
  Condition = ""  # fill this in as "Tumor" or "Normal"
)

# Save it as CSV
write_csv(sample_info, "sample_info_template.csv")
