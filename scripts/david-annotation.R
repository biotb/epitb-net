# Read input data

df_input <- read.csv("outliers_interactions.csv", header = TRUE, sep = ',')
df_outliers <- df_input[df_input$Gene_1_ENTREZid != "Not found in GFF" & df_input$Gene_2_ENTREZid != "Not found in GFF", ]

# Read DAVID enrichment 
df_david <- read.csv("outliers_DAVID_enrichment.csv", header = TRUE, sep = ',')
# Add GO columns

#Category	Term
#GOTERM_BP_DIRECT	GO:0009405~pathogenesis
#GOTERM_BP_DIRECT	GO:0071555~cell wall organization
#GOTERM_CC_DIRECT	GO:0005618~cell wall
#GOTERM_CC_DIRECT	GO:0005886~plasma membrane
#GOTERM_CC_DIRECT	GO:0005829~cytosol
#GOTERM_CC_DIRECT	GO:0005887~integral component of plasma membrane
#GOTERM_MF_DIRECT	GO:0005524~ATP binding
#GOTERM_MF_DIRECT	GO:0004721~phosphoprotein phosphatase activity

color_neutro = "#FFFFFF"

df_outliers['pathogenesis_1'] <- "NO"
df_outliers['pathogenesis_COLOR_1'] <- color_neutro 
df_outliers['cell_wall_organization_1'] <- "NO"
df_outliers['cell_wall_organization_COLOR_1'] <- color_neutro 
df_outliers['cell_wall_1'] <- "NO"
df_outliers['cell_wall_COLOR_1'] <- color_neutro 
df_outliers['plasma_membrane_1'] <- "NO"
df_outliers['plasma_membrane_COLOR_1'] <- color_neutro 
df_outliers['cytosol_1'] <- "NO"
df_outliers['cytosol_COLOR_1'] <- color_neutro 
df_outliers['integral_component_1'] <- "NO"
df_outliers['integral_component_COLOR_1'] <- color_neutro
df_outliers['ATP_binding_1'] <- "NO"
df_outliers['ATP_binding_COLOR_1'] <- color_neutro 
df_outliers['phosphoprotein_1'] <- "NO"
df_outliers['phosphoprotein_COLOR_1'] <- color_neutro

df_outliers['pathogenesis_2'] <- "NO"
df_outliers['pathogenesis_COLOR_2'] <- color_neutro 
df_outliers['cell_wall_organization_2'] <- "NO"
df_outliers['cell_wall_organization_COLOR_2'] <- color_neutro 
df_outliers['cell_wall_2'] <- "NO"
df_outliers['cell_wall_COLOR_2'] <- color_neutro 
df_outliers['plasma_membrane_2'] <- "NO"
df_outliers['plasma_membrane_COLOR_2'] <- color_neutro 
df_outliers['cytosol_2'] <- "NO"
df_outliers['cytosol_COLOR_2'] <- color_neutro 
df_outliers['integral_component_2'] <- "NO"
df_outliers['integral_component_COLOR_2'] <- color_neutro
df_outliers['ATP_binding_2'] <- "NO"
df_outliers['ATP_binding_COLOR_2'] <- color_neutro 
df_outliers['phosphoprotein_2'] <- "NO"
df_outliers['phosphoprotein_COLOR_2'] <- color_neutro

# Annotate table with GOs
gos = c("GO:0009405~pathogenesis", "GO:0071555~cell wall organization", "GO:0005618~cell wall", "GO:0005886~plasma membrane", "GO:0005829~cytosol", "GO:0005887~integral component of plasma membrane", "GO:0005524~ATP binding", "GO:0004721~phosphoprotein phosphatase activity")
cols = c("pathogenesis", "cell_wall_organization", "cell_wall", "plasma_membrane", "cytosol", "integral_component", "ATP_binding", "phosphoprotein")
cols_colors = c("pathogenesis_COLOR", "cell_wall_organization_COLOR", "cell_wall_COLOR", "plasma_membrane_COLOR", "cytosol_COLOR", "integral_component_COLOR", "ATP_binding_COLOR", "phosphoprotein_COLOR")
colors = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

for (i in 1:length(gos)) {
  genes <- df_david[df_david$Term == gos[i], ]$Genes
  # print(genes)
  list_genes <- strsplit(genes, ",", fixed=T)
  # print(list_genes)
  for (gene in list_genes[[1]]) {
    g <- gsub(" ","",gene)
    # print(g)
    col = cols[i]
    col_col = cols_colors[i]
    color = colors[i]
    # print(color)
    if(g %in% df_outliers$Gene_1_ENTREZid){
      col = paste(col, "_1", sep = "")
      # print(col)
      col_col = paste(col_col, "_1", sep = "")
      # print(col_col)
      df_outliers[grepl(g, df_outliers$Gene_1_ENTREZid), ][[col]] <- "YES"
      df_outliers[grepl(g, df_outliers$Gene_1_ENTREZid), ][[col_col]] <- color
    }
    col = cols[i]
    col_col = cols_colors[i]
    if(g %in% df_outliers$Gene_2_ENTREZid){
      col = paste(col, "_2", sep = "")
      # print(col)
      col_col = paste(col_col, "_2", sep = "")
      # print(col_col)
      df_outliers[grepl(g, df_outliers$Gene_2_ENTREZid), ][[col]] <- "YES"
      df_outliers[grepl(g, df_outliers$Gene_2_ENTREZid), ][[col_col]] <- color
    }
  }
  
}

## Save outliers with GO annotations
output <- df_outliers
write.csv(output, "outliers-gos.csv", row.names=FALSE, quote=FALSE)

## Generate plot
library(ggplot2)

df_david_plot <- df_david
df_david_plot$Term_plot <- paste(toupper(substr(gsub("^GO:[^~]+~", "", df_david_plot$Term), 1, 1)), substr(gsub("^GO:[^~]+~", "", df_david_plot$Term), 2, nchar(gsub("^GO:[^~]+~", "", df_david_plot$Term))), sep="")
names(df_david_plot)[names(df_david_plot)=="X."] <- "Percentage"

ggplot(df_david_plot, aes(reorder(Term_plot,
                                  Percentage),
                          Percentage, fill=Category)) +
  # 1) 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("#FED439B2", "#709AE1B2", "#D2AF81B2"), labels = c("BP", "CC", "MF")) +
  theme_light() +
  coord_flip() +
  geom_text(aes(label=sprintf(Pvalue, fmt = '%#.3f')), vjust=-0.3, size=2.5) +
  xlab("Enriched term") +
  ylab("% Genes")

ggsave("david-enrichment-outliers.png")
