# Create table with annotations for phandango plot with allele distribution.
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Run: Rscript phandango.R outliers.csv glgB
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Two arguments must be supplied: Input file, Gene name.n", call.=FALSE)
}

input_full_filepath <- args[1]
print(paste("Input file: ", input_full_filepath, sep = ""))

input_gene = args[2]
print(paste("Gene: ", input_gene, sep = ""))

# Read input data
print("Reading input data...")
df_input <- read.csv(input_full_filepath, header = TRUE, stringsAsFactors=FALSE)
print(dim(df_input))
print("   Done!")

# # Filter records by gene name and obtain dataframe for nucleotide extraction
library(dplyr)
print("Filtering genes...")
print("  Filtering gene position 1...")
df_filtered_1 <- df_input[(df_input$Gene_1_name == input_gene & df_input$Gene_2_name != "group"), ]
print(dim(df_filtered_1))
print("  Creating dataframe position 1...")
df_gene_1 <- df_filtered_1[, c("Position_1", "Gene_1_name", "Position_2", "Gene_2_name", "MI")]
colnames(df_gene_1) <- c("Position_1", "Gene_1_name", "Position_2", "Gene_2_name", "MI")
print(dim(df_gene_1))

print("  Filtering gene position 2...")
df_filtered_2 <- df_input[(df_input$Gene_2_name == input_gene & df_input$Gene_1_name != "group"), ]
print(dim(df_filtered_2))
print("  Creating dataframe position 2...")
df_gene_2 <- df_filtered_2[, c("Position_2", "Gene_2_name", "Position_1", "Gene_1_name", "MI")]
colnames(df_gene_2) <- c("Position_1", "Gene_1_name", "Position_2", "Gene_2_name", "MI")
print(dim(df_gene_2))

# Union gene dataframes
print("Union gene dataframes...")
df_gene <- rbind(df_gene_1, df_gene_2)
df_gene <- df_gene[order(df_gene$Position_1, -df_gene$MI, df_gene$Position_2), ]
print(dim(df_gene))

# Read fasta alignment with 254 strains with seqinr
# install.packages("seqinr")
# https://davetang.org/muse/2013/05/09/using-the-r-seqinr-package/
# https://stackoverflow.com/questions/45281606/parsing-information-out-of-sequencing-data
library("seqinr")
print("Reading fasta alignment...")
fastaFile <- read.fasta(file = "../concatenated.out")
print("   Done!")

# Read original phandango annotation file
df_phandango <- read.csv("phandango-annotations.csv", header = TRUE, stringsAsFactors=FALSE)
names(df_phandango)[names(df_phandango)=="Submitter.colour"] <- "Submitter:colour"

# Create phandango annotation file
print("Getting nucleotides from positions...")
pos_prev <- ""
for(row in 1:nrow(df_gene)){
  pos1 <- df_gene$Position_1[row]
  gene1 <- df_gene$Gene_1_name[row]
  if(pos1 == pos_prev){
    pos2 <- df_gene$Position_2[row]
    gene2 <- df_gene$Gene_2_name[row]
    column_name <- paste(gene2, "_", pos2, ":o1", sep = "")
    df_temp <- data.frame(df_phandango$name)
    colnames(df_temp) <- c("name")
    df_temp[[column_name]] <- NA
    pos_fasta <- getFrag(fastaFile, begin = pos2, end = pos2)
    for(i in 1:length(pos_fasta)){
      df_temp[df_temp$name == attr(pos_fasta[[i]],"seqMother"), column_name] <- ifelse(toupper(pos_fasta[[i]][1]) == "?", "-", toupper(pos_fasta[[i]][1]))
    }
    df_phandango <- cbind(df_phandango, df_temp[c(column_name)])
    
  } else {
    column_name <- paste(gene1, "_", pos1, ":o1", sep = "")
    df_temp <- data.frame(df_phandango$name)
    colnames(df_temp) <- c("name")
    df_temp[[column_name]] <- NA
    pos_fasta <- getFrag(fastaFile, begin = pos1, end = pos1)
    for(i in 1:length(pos_fasta)){
      df_temp[df_temp$name == attr(pos_fasta[[i]],"seqMother"), column_name] <- ifelse(toupper(pos_fasta[[i]][1]) == "?", "-", toupper(pos_fasta[[i]][1]))
    }
    df_phandango <- cbind(df_phandango, df_temp[c(column_name)])

    pos2 <- df_gene$Position_2[row]
    gene2 <- df_gene$Gene_2_name[row]
    column_name <- paste(gene2, "_", pos2, ":o1", sep = "")
    df_temp <- data.frame(df_phandango$name)
    colnames(df_temp) <- c("name")
    df_temp[[column_name]] <- NA
    pos_fasta <- getFrag(fastaFile, begin = pos2, end = pos2)
    for(i in 1:length(pos_fasta)){
      df_temp[df_temp$name == attr(pos_fasta[[i]],"seqMother"), column_name] <- ifelse(toupper(pos_fasta[[i]][1]) == "?", "-", toupper(pos_fasta[[i]][1]))
    }
    df_phandango <- cbind(df_phandango, df_temp[c(column_name)])
    pos_prev <- pos1
  }
print(paste(" Row:", row))
}
print("   Done!")

print("Saving phandango annotation file...")
# Save new phangando annotation file
write.csv(df_phandango, file = paste("phandango-annotations-", input_gene, ".csv", sep = ""), row.names = FALSE)

# Add MDR and XDR to phandango annotation file
## Read phandango annotation file
df_phandango <- read.csv(paste("phandango-annotations-", input_gene, ".csv", sep = ""), header = TRUE, stringsAsFactors=FALSE, check.names = FALSE)
names(df_phandango)[names(df_phandango)=="Submitter.colour"] <- "Submitter:colour"
for (i in 4:length(colnames(df_phandango))) {
  colnames(df_phandango)[i] <- sub("\\.", ":", colnames(df_phandango)[i])
}

## Read MDR and XDR iTOL annotations
df_itol <- read.csv("dataset_binary_DATA_MDR_XDR.csv", header = TRUE, stringsAsFactors=FALSE)

## Add GAR (MDR or XDR) annotations to phandango annotation file
df_phandango_final <- df_phandango
column_color <- "GAR:colour"
df_phandango_final$GAR <- ""
df_phandango_final[[column_color]] <- "#D9D9D9"
for (g in df_phandango_final$name) {
  if(df_itol[df_itol$genome == g, "MDR"] == "1"){
    df_phandango_final[df_phandango_final$name == g, "GAR"] <- "MDR"
    df_phandango_final[df_phandango_final$name == g, "GAR:colour"] <- "#93c47d"
  }
  if(df_itol[df_itol$genome == g, "XDR"] == "1"){
    df_phandango_final[df_phandango_final$name == g, "GAR"] <- "XDR"
    df_phandango_final[df_phandango_final$name == g, "GAR:colour"] <- "#e06666"
  }
}
df_phandango_final <- df_phandango_final[, c("name", "Submitter", "Submitter:colour", "GAR", column_color)]
for (i in 4:length(colnames(df_phandango))) {
  df_phandango_final <- cbind(df_phandango_final, df_phandango[i])
}

# Save new phangando annotation file
write.csv(df_phandango_final, file = paste("phandango-annotations-MDR-XDR-", input_gene, ".csv", sep = ""), row.names = FALSE)



