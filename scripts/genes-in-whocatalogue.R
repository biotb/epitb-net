# Look for genes in WHO catalogue and annotation with WHO drugs.

path = "/home/"

# Read outliers with GO annotations
df_outliers_gos <- read.csv(paste(path, "outliers-gos.csv", sep = ""), header = TRUE, sep = ',')

# Create dataframes for genes
## Read tables with genes connected to oppA and glgB exported from Cytoscape
df_oppa <- read.csv(paste(path, "oppA.csv", sep = ""), header = TRUE, sep = ',')
df_glgb <- read.csv(paste(path, "glgB.csv", sep = ""), header = TRUE, sep = ',')

# Annotate genes of glgB and oppA with product, GO, WHO-drug and WHO-position 
## Read table with gene and product
df_product <- read.csv(paste(path, "gene_product.csv", sep = ""), header = TRUE, sep = ',')

## Read WHO catalogue
df_who_catalogue <- read.csv(paste(path, "WHO-UCN-GTB-PCI-2021.7-eng.csv", sep = ""), header = TRUE, sep = ',')

## Create supplementary table of genes
### Create function
get_gene_table = function(df_g, df_p, df_go) {
  df_gw <- data.frame(df_g$name)
  print(dim(df_gw))
  names(df_gw)[names(df_gw)=="df_g.name"] <- "gene"
  df_gw$pathogenesis <- ""
  df_gw$cell_wall_organization <- ""
  df_gw$cell_wall <- ""
  df_gw$plasma_membrane <- ""
  df_gw$cytosol <- ""
  df_gw$integral_component <- ""
  df_gw$ATP_binding <- ""
  df_gw$phosphoprotein <- ""
  df_gw$product <- ""
  df_gw$drugs <- ""
  df_gw$genome_position <- ""
  for (i in 1:nrow(df_gw)) {
    gene <- df_gw$gene[i]
    if(gene %in% df_p$gene){
      df_gw$product[i] <- df_p[df_p$gene == gene, ]$product
    } else {
      df_gw$product[i] <- "pendiente"
    }
    df_gwene_who <- df_who_catalogue[grepl(gene, df_who_catalogue[["variant..common_name."]]), c("drug", "variant..common_name.", "Genome.position")]
    if(nrow(df_gwene_who) > 0){
      drugs <- unique(df_gwene_who$drug)
      genome_position <- unique(df_gwene_who$Genome.position)
      df_gw[df_gw$gene == gene, ]$drugs <- paste(drugs[!is.na(drugs)], collapse = ",")
      df_gw[df_gw$gene == gene, ]$genome_position <- paste(genome_position[!is.na(genome_position)], collapse = ",")
    }
  }
  return(df_gw)
}

gene = "glgB"
df_genes_glgb <- get_gene_table(df_glgb, df_product, df_outliers_gos)
write.csv(df_genes_glgb, file = paste(path, "who-drug-", gene, ".csv", sep = ""), row.names = FALSE)
df_drug_glgb <- df_genes_glgb[df_genes_glgb$drugs != "", c("gene", "drugs")]
write.csv(df_drug_glgb, file = paste(path, "drug-", gene, ".csv", sep = ""), row.names = FALSE)

gene = "oppA"
df_genes_oppa <- get_gene_table(df_oppa, df_product, df_outliers_gos)
write.csv(df_genes_oppa, file = paste(path, "who-drug-", gene , ".csv", sep = ""), row.names = FALSE)
df_drug_oppa <- df_genes_oppa[df_genes_oppa$drugs != "", c("gene", "drugs")]
write.csv(df_drug_oppa, file = paste(path, "drug-", gene, ".csv", sep = ""), row.names = FALSE)
