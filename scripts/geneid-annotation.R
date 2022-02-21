# Gene ENTREZ id annotation with GFF file

path = "/home/"

# Read file with genename annotation
input <- read.csv(paste(path, "genenames.csv", sep = ""), header = TRUE, sep = ',')
print(paste("Total de asociaciones outliers: ", dim(input)[1]))

# Filter genes without groups (group_)
input_outliers <- input[!grepl("group", input$Gene_1) & !grepl("group", input$Gene_2), ]
print(paste("Total de asociaciones outliers sin groups: ", dim(input_outliers)[1]))

# Get unique genes
library(dplyr)

get_unique_genes <- function(df_data) {
  df_gene1 <- data.frame(df_data[, c("Gene_1")])
  colnames(df_gene1) <- c("Gene")
  df_gene2 <- data.frame(df_data[, c("Gene_2")])
  colnames(df_gene2) <- c("Gene")
  df_genes <- rbind(df_gene1, df_gene2)

  df_genes <- data.frame(df_genes[!duplicated(df_genes[ , c("Gene")]), ])
  colnames(df_genes) <- c("Gene")
  df_genes <- data.frame(df_genes[order(df_genes$Gene),])
  colnames(df_genes) <- c("Gene")
  return(df_genes)
}

df_unique_genes_outliers <- get_unique_genes(input_outliers)

# Obtain Gene ID from gff file
# BiocManager::install("Biostrings")
library("biomartr")
library(stringr)

# Load GFF file
gff_file <- biomartr::read_gff(file = paste(path, "GCF_000195955.2_ASM19595v2_genomic.gff", sep = ""))

# Filter genes
df_gff_file <- gff_file[gff_file$type == "gene", ]

# Extract gene_name and gene_id (GeneID) as new columns
# ID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;
for (i in 1:nrow(df_gff_file)) {
  df_gff_file$gene_name[i] <- str_match(df_gff_file$attribute[i],"Name=([^;]+);")[1,2]
  df_gff_file$gene_id[i] <- str_match(df_gff_file$attribute[i],"GeneID:([^;]+);")[1,2]
}

# Annotate genes with Gene ID using new columns from gff
# We removed the _N suffix from the gene names
not_found = 0
found = 0
for(i in 1:nrow(df_unique_genes_outliers)) {
  if(df_unique_genes_outliers$Gene[i] %in% df_gff_file$gene_name) {
    found = found + 1
    df_unique_genes_outliers$gene_id[i] = df_gff_file$gene_id[df_gff_file$gene_name == df_unique_genes_outliers$Gene[i]]
  } else {
    not_found = not_found + 1
    df_unique_genes_outliers$gene_id[i] = "Not found in GFF"
  }
}
print("df_unique_genes_outliers")
print(paste("Genes found: ", found))
print(paste("Genes not found: ", not_found))

# Outliers annotation with gene id
input_outliers$Gene_1_id <- "Not found in GFF"
input_outliers$Gene_2_id <- "Not found in GFF"

for(i in 1:nrow(df_unique_genes_outliers)) {
  gene = df_unique_genes_outliers$Gene[i]
  gene_id = df_unique_genes_outliers$gene_id[i]
  if(gene %in% input_outliers$Gene_1){
    input_outliers[input_outliers$Gene_1 == gene, ]$Gene_1_id <- gene_id
  }
  if(gene %in% input_outliers$Gene_2){
    input_outliers[input_outliers$Gene_2 == gene, ]$Gene_2_id <- gene_id
  }
}
output <- input_outliers[, c("pos1", "Gene_1_id", "gene1", "pos2", "Gene_2_id", "gene2", "distance", "type", "mi", "extreme_outlier")]
colnames(output) <- c("Position_1", "Gene_1_ENTREZid", "Gene_1_name", "Position_2", "Gene_2_ENTREZid", "Gene_2_name", "Distance", "Type", "MI", "Extreme_outlier")
write.csv(output, paste(path, "geneids.csv", sep = ""), row.names=FALSE, quote=FALSE)
