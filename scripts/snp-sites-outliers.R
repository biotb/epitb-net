# Compare positions from outliers to SNP-sites.

# Load outlier interactions
input_full_filepath <- "/home/outliers-gos.csv"
input <- read.csv(input_full_filepath)

# Create list of unique positions (Position_1 and Position_2)
df_positions <- as.vector(as.matrix(input[, c("Position_1", "Position_2")]))
vec_unique_positions <- sort(unique(df_positions))

# Open VCF file creted by SNP-sites
library(vcfR)
vcf <- read.vcfR("/home/amas-roary.vcf", verbose = FALSE)
all_pos = getPOS(vcf)
positions_found <- integer()
print(paste("Total SpydrPick unique positions:", length(vec_unique_positions)))
for (i in 1:length(vec_unique_positions)) {
  pos <- vec_unique_positions[i]
  if(pos %in% all_pos){
    # print(paste("Position found!", pos))
    positions_found <- c(positions_found, pos)
  }
}

print(paste("Total unique positions from outlier interactions found in SNPs from SNP-sites:", length(positions_found)))
