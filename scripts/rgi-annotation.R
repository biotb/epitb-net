# Strain annotation with genotypic antibiotic resistance.

## Read 254 genome list
input_full_filepath <- "/home/rgi_card/genome-list-complete.txt"
df_genomes <- read.csv(input_full_filepath, header = FALSE, stringsAsFactors=FALSE)
colnames(df_genomes) <- c("genome")

## Read Resistance Gene Identifier txt files from RIG-CARD output
### Get antibiotics of interest
# MDR: isoniazid and rifampicin
# XDR: isoniazid and rifampicin and fluoroquinolone and (kanamycin or amikacin or capreomycin)
df_genomes_rig <- data.frame(genome=character(),
                             isoniazid=character(),
                             rifampicin=character(),
                             fluoroquinolone=character(),
                             kanamycin=character(),
                             amikacin=character(),
                             capreomycin=character(),
                stringsAsFactors=FALSE)
input_path <- "/home/rgi_card/"
for (g in df_genomes$genome) {
  input_file <- paste(g, "_genomic.fna.rig.txt", sep = "")
  input_full_filepath <- paste(input_path, input_file, sep = "")
  print(paste("Processing:", g))
  df_rig_file <- read.csv(input_full_filepath, header = TRUE, stringsAsFactors=FALSE, sep = "\t")
  new_row <- data.frame(genome = g, isoniazid = "FALSE",
                             rifampicin = "FALSE",
                             fluoroquinolone = "FALSE",
                             kanamycin = "FALSE",
                             amikacin = "FALSE",
                             capreomycin = "FALSE")
  df_genomes_rig <- rbind(df_genomes_rig, new_row)
  if(length(grep("isoniazid", df_rig_file[["Best_Hit_ARO"]])) > 0){
    df_genomes_rig[df_genomes_rig$genome == g, ]$isoniazid = "TRUE"
  }
  if(length(grep("rifampicin", df_rig_file[["Best_Hit_ARO"]])) > 0){
    df_genomes_rig[df_genomes_rig$genome == g, ]$rifampicin = "TRUE"
  }
  if(length(grep("fluoroquinolone", df_rig_file[["Best_Hit_ARO"]])) > 0){
    df_genomes_rig[df_genomes_rig$genome == g, ]$fluoroquinolone = "TRUE"
  }
  if(length(grep("kanamycin", df_rig_file[["Best_Hit_ARO"]])) > 0){
    df_genomes_rig[df_genomes_rig$genome == g, ]$kanamycin = "TRUE"
  }
  if(length(grep("amikacin", df_rig_file[["Best_Hit_ARO"]])) > 0){
    df_genomes_rig[df_genomes_rig$genome == g, ]$amikacin = "TRUE"
  }
  if(length(grep("capreomycin", df_rig_file[["Best_Hit_ARO"]])) > 0){
    df_genomes_rig[df_genomes_rig$genome == g, ]$capreomycin = "TRUE"
  }
}

## Create dataset_binary DATA for iTOL
df_data_itol <- data.frame(genome=character(), MDR=character(), XDR=character(),
                stringsAsFactors=FALSE)
for (g in df_genomes_rig$genome) {
  gs <- sub("\\..+$", "", g)
  print(gs)
  new_row <- data.frame(genome=gs, MDR="-1", XDR="-1")
  df_data_itol <- rbind(df_data_itol, new_row)
  if(nrow(df_genomes_rig[df_genomes_rig$genome == g & df_genomes_rig$isoniazid == "TRUE" & df_genomes_rig$rifampicin == "TRUE", ]) > 0){
    print("  MDR!!")
    df_data_itol[df_data_itol$genome == gs, ]$MDR = "1"
  }
  if(nrow(df_genomes_rig[df_genomes_rig$genome == g & df_genomes_rig$isoniazid == "TRUE" & df_genomes_rig$rifampicin == "TRUE" & df_genomes_rig$fluoroquinolone == "TRUE" & (df_genomes_rig$kanamycin == "TRUE" | df_genomes_rig$amikacin == "TRUE" | df_genomes_rig$capreomycin == "TRUE"), ]) > 0){
    print("  XDR!!")
    df_data_itol[df_data_itol$genome == gs, ]$XDR = "1"
  }
}

output_full_filepath <- "/home/dataset_binary_DATA_MDR_XDR.csv"
write.table(df_data_itol,output_full_filepath, row.names=FALSE, quote=FALSE, sep = ",")

## Percentage of MDR and XDR
df_data <- read.csv("/home/dataset_binary_DATA_MDR_XDR.csv", header = TRUE, stringsAsFactors=FALSE)
print((nrow(df_data[df_data$MDR == 1, ])/nrow(df_data))*100)
pint((nrow(df_data[df_data$XDR == 1, ])/nrow(df_data))*100)
