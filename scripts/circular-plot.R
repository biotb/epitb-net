# Read outliers with GO annotations

df_outliers <- read.csv("outliers-gos.csv", header = TRUE, sep = ',')

## Scaled mi between 0 and 1
# Scaling 0 to 1
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_outliers['scaled01_mi'] <- scale01(df_outliers['MI'])

# Create dataframes for genes
## Read gene tables exported from Cytoscape
df_oppa <- read.csv("oppA.csv", header = TRUE, sep = ',')
df_glgb <- read.csv("glgB.csv", header = TRUE, sep = ',')

## Create dataframes for each module for each GO
### BP
df_outliers_BP_oppa <- df_outliers[(df_outliers$pathogenesis_1 == "YES" | df_outliers$cell_wall_organization_1 == "YES") &  (df_outliers$pathogenesis_2 == "YES" | df_outliers$cell_wall_organization_2 == "YES") & (df_outliers$Gene_1_name %in% df_oppa$name | df_outliers$Gene_2_name %in% df_oppa$name), ]
df_outliers_BP_glgb <- df_outliers[(df_outliers$pathogenesis_1 == "YES" | df_outliers$cell_wall_organization_1 == "YES") &  (df_outliers$pathogenesis_2 == "YES" | df_outliers$cell_wall_organization_2 == "YES") & (df_outliers$Gene_1_name %in% df_glgb$name | df_outliers$Gene_2_name %in% df_glgb$name), ]

### BP adding genes 
df_outliers_BP_oppa_2 <- df_outliers[(df_outliers$Gene_1_name == "oppA" & df_outliers$Gene_2_name %in% df_outliers_BP_oppa$Gene_2_name) | (df_outliers$Gene_2_name == "oppA" & df_outliers$Gene_1_name %in% df_outliers_BP_oppa$Gene_1_name), ]
df_outliers_BP_oppa_all <- rbind(df_outliers_BP_oppa, df_outliers_BP_oppa_2)

df_outliers_BP_glgb_2 <- df_outliers[(df_outliers$Gene_1_name == "glgB" & df_outliers$Gene_2_name %in% df_outliers_BP_glgb$Gene_2_name) | (df_outliers$Gene_2_name == "glgB" & df_outliers$Gene_1_name %in% df_outliers_BP_glgb$Gene_1_name), ]
df_outliers_BP_glgb_all <- rbind(df_outliers_BP_glgb, df_outliers_BP_glgb_2)

### CC
df_outliers_CC_oppa <- df_outliers[(df_outliers$cell_wall_1 == "YES" | df_outliers$plasma_membrane_1 == "YES" | df_outliers$cytosol_1 == "YES" | df_outliers$integral_component_1 == "YES") & (df_outliers$cell_wall_2 == "YES" | df_outliers$plasma_membrane_2 == "YES" | df_outliers$cytosol_2 == "YES" | df_outliers$integral_component_2 == "YES") & (df_outliers$Gene_1_name %in% df_oppa$name | df_outliers$Gene_2_name %in% df_oppa$name), ]
df_outliers_CC_glgb <- df_outliers[(df_outliers$cell_wall_1 == "YES" | df_outliers$plasma_membrane_1 == "YES" | df_outliers$cytosol_1 == "YES" | df_outliers$integral_component_1 == "YES") & (df_outliers$cell_wall_2 == "YES" | df_outliers$plasma_membrane_2 == "YES" | df_outliers$cytosol_2 == "YES" | df_outliers$integral_component_2 == "YES") & (df_outliers$Gene_1_name %in% df_glgb$name | df_outliers$Gene_2_name %in% df_glgb$name), ]

### CC adding genes 
df_outliers_CC_oppa_2 <- df_outliers[(df_outliers$Gene_1_name == "oppA" & df_outliers$Gene_2_name %in% df_outliers_CC_oppa$Gene_2_name) | (df_outliers$Gene_2_name == "oppA" & df_outliers$Gene_1_name %in% df_outliers_CC_oppa$Gene_1_name), ]
df_outliers_CC_oppa_all <- rbind(df_outliers_CC_oppa, df_outliers_CC_oppa_2)

df_outliers_CC_glgb_2 <- df_outliers[(df_outliers$Gene_1_name == "glgB" & df_outliers$Gene_2_name %in% df_outliers_CC_glgb$Gene_2_name) | (df_outliers$Gene_2_name == "glgB" & df_outliers$Gene_1_name %in% df_outliers_CC_glgb$Gene_1_name), ]
df_outliers_CC_glgb_all <- rbind(df_outliers_CC_glgb, df_outliers_CC_glgb_2)

### MF
df_outliers_MF_oppa <- df_outliers[(df_outliers$ATP_binding_1 == "YES" | df_outliers$phosphoprotein_1 == "YES") & (df_outliers$ATP_binding_2 == "YES" | df_outliers$phosphoprotein_2 == "YES") & (df_outliers$Gene_1_name %in% df_oppa$name | df_outliers$Gene_2_name %in% df_oppa$name), ]
df_outliers_MF_glgb <- df_outliers[(df_outliers$ATP_binding_1 == "YES" | df_outliers$phosphoprotein_1 == "YES") & (df_outliers$ATP_binding_2 == "YES" | df_outliers$phosphoprotein_2 == "YES") & (df_outliers$Gene_1_name %in% df_glgb$name | df_outliers$Gene_2_name %in% df_glgb$name), ]

### MF adding genes 
df_outliers_MF_oppa_2 <- df_outliers[(df_outliers$Gene_1_name == "oppA" & df_outliers$Gene_2_name %in% df_outliers_MF_oppa$Gene_2_name) | (df_outliers$Gene_2_name == "oppA" & df_outliers$Gene_1_name %in% df_outliers_MF_oppa$Gene_1_name), ]
df_outliers_MF_oppa_all <- rbind(df_outliers_MF_oppa, df_outliers_MF_oppa_2)

df_outliers_MF_glgb_2 <- df_outliers[(df_outliers$Gene_1_name == "glgB" & df_outliers$Gene_2_name %in% df_outliers_MF_glgb$Gene_2_name) | (df_outliers$Gene_2_name == "glgB" & df_outliers$Gene_1_name %in% df_outliers_MF_glgb$Gene_1_name), ]
df_outliers_MF_glgb_all <- rbind(df_outliers_MF_glgb, df_outliers_MF_glgb_2)

# Create function to get unique genes with GOs and colors and positions

get_genes_unique = function(df_out) {
  dfgene1 <- df_out[, c("Gene_1_name", "pathogenesis_1", "pathogenesis_COLOR_1", "cell_wall_organization_1", "cell_wall_organization_COLOR_1", "cell_wall_1", "cell_wall_COLOR_1", "plasma_membrane_1", "plasma_membrane_COLOR_1", "cytosol_1", "cytosol_COLOR_1", "integral_component_1", "integral_component_COLOR_1", "ATP_binding_1", "ATP_binding_COLOR_1", "phosphoprotein_1", "phosphoprotein_COLOR_1")]
  colnames(dfgene1) <- c("gene", "pathogenesis", "pathogenesis_COLOR", "cell_wall_organization", "cell_wall_organization_COLOR", "cell_wall", "cell_wall_COLOR", "plasma_membrane", "plasma_membrane_COLOR", "cytosol", "cytosol_COLOR", "integral_component", "integral_component_COLOR", "ATP_binding", "ATP_binding_COLOR", "phosphoprotein", "phosphoprotein_COLOR")

  dfgene2 <- df_out[, c("Gene_2_name", "pathogenesis_2", "pathogenesis_COLOR_2", "cell_wall_organization_2", "cell_wall_organization_COLOR_2", "cell_wall_2", "cell_wall_COLOR_2", "plasma_membrane_2", "plasma_membrane_COLOR_2", "cytosol_2", "cytosol_COLOR_2", "integral_component_2", "integral_component_COLOR_2", "ATP_binding_2", "ATP_binding_COLOR_2", "phosphoprotein_2", "phosphoprotein_COLOR_2")]
  colnames(dfgene2) <- c("gene", "pathogenesis", "pathogenesis_COLOR", "cell_wall_organization", "cell_wall_organization_COLOR", "cell_wall", "cell_wall_COLOR", "plasma_membrane", "plasma_membrane_COLOR", "cytosol", "cytosol_COLOR", "integral_component", "integral_component_COLOR", "ATP_binding", "ATP_binding_COLOR", "phosphoprotein", "phosphoprotein_COLOR")

  df_genes_unique <- rbind(dfgene1, dfgene2)
  df_genes_unique <- df_genes_unique[,c("gene", "pathogenesis", "pathogenesis_COLOR", "cell_wall_organization", "cell_wall_organization_COLOR", "cell_wall", "cell_wall_COLOR", "plasma_membrane", "plasma_membrane_COLOR", "cytosol", "cytosol_COLOR", "integral_component", "integral_component_COLOR", "ATP_binding", "ATP_binding_COLOR", "phosphoprotein", "phosphoprotein_COLOR")]
  df_genes_unique <- df_genes_unique %>% distinct()
  
  df_genes_unique['start'] <- 0
  df_genes_unique['end'] <- 0
  
  # Load AMAS partitions.txt file
  partitions <- read.csv("partitions.txt", header = FALSE, sep = '=')
  colnames(partitions) <- c("gene","range")
  df_genes <- partitions %>%
    separate(gene, c("p", "gene_suffix"), "_", remove = FALSE, extra = "merge")
  
  df_genes$gene <- sub("_[0-9]+","",df_genes$gene_suffix)
  
  df_genes <- df_genes %>%
    separate(range, c("start", "end"), "-")
  
  df_genes$start <- as.numeric(as.character(df_genes$start))
  df_genes$end <- as.numeric(as.character(df_genes$end))
  
  df_genes$gene <- trimws(df_genes$gene, which = c("right"))

# Get start and end gene position from partitions.txt file
  for(i in 1:nrow(df_genes)) {
        gene <- df_genes$gene[i]
        df_range <- df_genes[df_genes$gene == gene, ]
        start = min(df_genes[df_genes$gene == gene, c("start")])
        end = max(df_genes[df_genes$gene == gene, c("end")])
        if(gene %in% df_genes_unique$gene) {
          df_genes_unique[df_genes_unique$gene == gene, ]$start <- start
          df_genes_unique[df_genes_unique$gene == gene, ]$end <- end
        }
  }
  # https://statisticsglobe.com/convert-data-frame-column-to-numeric-in-r
  df_genes_unique$start <- as.double(df_genes_unique$start)
  df_genes_unique$end <- as.double(df_genes_unique$end)
  df_genes_unique <- df_genes_unique[order(df_genes_unique$start),]
  
  return(df_genes_unique)
}

## Create unique genes for each GO and module

library(dplyr)
library(tidyr)

### BP
df_genes_unique_BP_oppa <- get_genes_unique(df_outliers_BP_oppa_all)
df_genes_unique_BP_glgb <- get_genes_unique(df_outliers_BP_glgb_all)

### CC
df_genes_unique_CC_oppa <- get_genes_unique(df_outliers_CC_oppa_all)
df_genes_unique_CC_glgb <- get_genes_unique(df_outliers_CC_glgb_all)

### MF
df_genes_unique_MF_oppa <- get_genes_unique(df_outliers_MF_oppa_all)
df_genes_unique_MF_glgb <- get_genes_unique(df_outliers_MF_glgb_all)

## Define functions for circular plot With SEPARATED GOs
### Create functions

library(circlize)

`%notin%` <- Negate(`%in%`)
mis <- unique(df_outliers["scaled01_mi"])
colorred = "#E06666"

circlize_plot_BP = function(df_interactions, df_gu, gene) {
  factors <- df_gu[, "gene"]
  factors_col <- ifelse(factors == gene, "red", "black")
  xlim = data.matrix(df_gu[, c("start", "end")])
  # Size of gene names
  vcex = 0.5
  circos.clear()
  circos.par(cell.padding = c(0.02, 0, 0.02, 0), start.degree = 90, gap.degree = 0.5)
  circos.initialize(factors = factors, xlim = xlim)
  circos.track(factors = factors, ylim = c(0, 1), track.height = 0.05, track.margin=c(0,0.3), bg.col = df_gu$pathogenesis_COLOR,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(6, "mm"), 
                             CELL_META$sector.index, facing = "clockwise", cex = vcex, niceFacing = TRUE, col = factors_col[CELL_META$sector.numeric.index])
              }, bg.border=T, bg.lwd = .3)
  set_track_gap(mm_h(0.2)) # 2mm
  # Add new track for cell_wall_organization
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = df_gu$cell_wall_organization_COLOR, bg.border=T, bg.lwd = .3)

  for (row in 1:nrow(df_interactions)) {
    g1 <- df_interactions[row, "Gene_1_name"]
    p1 <- df_interactions[row, "Position_1"]
    g2 <- df_interactions[row, "Gene_2_name"]
    p2 <- df_interactions[row, "Position_2"]
    mi <- df_interactions[row, "scaled01_mi"]
    if((g1 == gene) | (g2 == gene)){
      c <- colorred
    } else{
      c <- "darkgray"
    }
    circos.link(g1, p1, g2, p2, lwd = .3, col = c, directional = 0)
  }
}

circlize_plot_CC = function(df_interactions, df_gu, gene) {
  factors <- df_gu[, "gene"]
  factors_col <- ifelse(factors == gene, "red", "black")
  xlim = data.matrix(df_gu[, c("start", "end")])
  # Size of gene names
  vcex = 0.5
  circos.clear()
  circos.par(cell.padding = c(0.02, 0, 0.02, 0), start.degree = 90, gap.degree = 0.5)
  circos.initialize(factors = factors, xlim = xlim)
  circos.track(factors = factors, ylim = c(0, 1), track.height = 0.05, track.margin=c(0,0.3), bg.col = df_gu$cell_wall_COLOR,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(6, "mm"), 
                             CELL_META$sector.index, facing = "clockwise", cex = vcex, niceFacing = TRUE, col = factors_col[CELL_META$sector.numeric.index])
              }, bg.border=T, bg.lwd = .3)
  set_track_gap(mm_h(0.2)) # 2mm
  # Add new track for cell_wall
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = df_gu$plasma_membrane_COLOR, bg.border=T, bg.lwd = .3)
  set_track_gap(mm_h(0.2)) # 2mm
  # Add new track for cytosol
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = df_gu$cytosol_COLOR, bg.border=T, bg.lwd = .3)
  set_track_gap(mm_h(0.2)) # 2mm
  # Add new track for integral_component
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = df_gu$integral_component_COLOR, bg.border=T, bg.lwd = .3)

  for (row in 1:nrow(df_interactions)) {
    g1 <- df_interactions[row, "Gene_1_name"]
    p1 <- df_interactions[row, "Position_1"]
    g2 <- df_interactions[row, "Gene_2_name"]
    p2 <- df_interactions[row, "Position_2"]
    mi <- df_interactions[row, "scaled01_mi"]
    if((g1 == gene) | (g2 == gene)){
      c <- colorred
    } else{
      c <- "darkgray"
    }
    circos.link(g1, p1, g2, p2, lwd = .3, col = c, directional = 0)
  }
}

circlize_plot_MF = function(df_interactions, df_gu, gene) {
  factors <- df_gu[, "gene"]
  factors_col <- ifelse(factors == gene, "red", "black")
  xlim = data.matrix(df_gu[, c("start", "end")])
  # Size of gene names
  vcex = 0.5
  circos.clear()
  circos.par(cell.padding = c(0.02, 0, 0.02, 0), start.degree = 90, gap.degree = 0.5)
  circos.initialize(factors = factors, xlim = xlim)
  circos.track(factors = factors, ylim = c(0, 1), track.height = 0.05, track.margin=c(0,0.3), bg.col = df_gu$ATP_binding_COLOR,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(6, "mm"), 
                             CELL_META$sector.index, facing = "clockwise", cex = vcex, niceFacing = TRUE, col = factors_col[CELL_META$sector.numeric.index])
              }, bg.border=T, bg.lwd = .3)
  set_track_gap(mm_h(0.2)) # 2mm
  # Add new track for phosphoprotein
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.col = df_gu$phosphoprotein_COLOR, bg.border=T, bg.lwd = .3)

  for (row in 1:nrow(df_interactions)) {
    g1 <- df_interactions[row, "Gene_1_name"]
    p1 <- df_interactions[row, "Position_1"]
    g2 <- df_interactions[row, "Gene_2_name"]
    p2 <- df_interactions[row, "Position_2"]
    mi <- df_interactions[row, "scaled01_mi"]
    #*c <- as.character(hash_colors[hash_colors$x == as.character(mi), "val"])
    if((g1 == gene) | (g2 == gene)){
      c <- colorred
    } else{
      c <- "darkgray"
    }
    circos.link(g1, p1, g2, p2, lwd = .3, col = c, directional = 0)
  }
}

## Define vertical legends for each GO and module

library(dplyr)
library(ComplexHeatmap)
# https://jokergoo.github.io/ComplexHeatmap-reference/book/
# sudo apt-get install libcairo2-dev
# BiocManager::install("ComplexHeatmap")

### Define functions to create legends

get_legend_BP = function(df_genes_unique) {
  lgd_points = Legend(at = c("Pathogenesis", "Cell wall organization"), type = "points", 
    legend_gp = gpar(col = c(unique(df_genes_unique[df_genes_unique$pathogenesis == "YES", ]$pathogenesis_COLOR), unique(df_genes_unique[df_genes_unique$cell_wall_organization == "YES", ]$cell_wall_organization_COLOR))), title_position = "topleft", 
    title = "Biological process", labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 7))
  lgd_list_vertical = packLegend(lgd_points)
  return(lgd_list_vertical)
}

get_legend_CC = function(df_genes_unique) {
  lgd_points = Legend(at = c("Cell wall", "Plasma membrane", "Cytosol", "Integral component plasma membr."), type = "points", 
    legend_gp = gpar(col = c(
      unique(df_genes_unique[df_genes_unique$cell_wall == "YES", ]$cell_wall_COLOR),
      unique(df_genes_unique[df_genes_unique$plasma_membrane == "YES", ]$plasma_membrane_COLOR),
      unique(df_genes_unique[df_genes_unique$cytosol == "YES", ]$cytosol_COLOR),
      unique(df_genes_unique[df_genes_unique$integral_component == "YES", ]$integral_component_COLOR)
      )), title_position = "topleft", 
    title = "Cellular component", labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 7))
  lgd_list_vertical = packLegend(lgd_points)
  return(lgd_list_vertical)
}

get_legend_MF = function(df_genes_unique) {
  lgd_points = Legend(at = c("ATP binding", "Phosphoprotein phosphatase activ."), type = "points", 
    legend_gp = gpar(col = c(
      unique(df_genes_unique[df_genes_unique$ATP_binding == "YES", ]$ATP_binding_COLOR),
      unique(df_genes_unique[df_genes_unique$phosphoprotein == "YES", ]$phosphoprotein_COLOR)
      )), title_position = "topleft", 
    title = "Molecular function", labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 7))
  lgd_list_vertical = packLegend(lgd_points)
  return(lgd_list_vertical)
}

### Create legend for each gene and GO
#### BP
#### oppA
lgd_list_vertical_BP_oppa <- get_legend_BP(df_genes_unique_BP_oppa)
#### glgB
lgd_list_vertical_BP_glgb <- get_legend_BP(df_genes_unique_BP_glgb)

#### CC
#### oppA
lgd_list_vertical_CC_oppa <- get_legend_CC(df_genes_unique_CC_oppa)
#### glgB
lgd_list_vertical_CC_glgb <- get_legend_CC(df_genes_unique_CC_glgb)

#### MF
#### oppA
lgd_list_vertical_MF_oppa <- get_legend_MF(df_genes_unique_MF_oppa)
#### glgB
lgd_list_vertical_MF_glgb <- get_legend_MF(df_genes_unique_MF_glgb)

## Create circular plots for each gene and GO with vertical legends
### Create functions for each GO

create_plot_BP = function(path, df_out_BP, df_genes_uni_BP, lgd_BP, g) {
  png(file=output_full_filepath, width = 1800, height = 1400, res = 300)
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
      just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  circlize_plot_BP(df_out_BP, df_genes_uni_BP, g)
  upViewport()
  draw(lgd_BP, x = circle_size, just = c("left", "top"))
  dev.off()
}

create_plot_CC = function(path, df_out_CC, df_genes_uni_CC, lgd_CC, g) {
  png(file=output_full_filepath, width = 1800, height = 1400, res = 300)
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
      just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  circlize_plot_CC(df_out_CC, df_genes_uni_CC, g)
  upViewport()
  draw(lgd_CC, x = circle_size, just = c("left", "top"))
  dev.off()
}

create_plot_MF = function(path, df_out_MF, df_genes_uni_MF, lgd_MF, g) {
  png(file=output_full_filepath, width = 1800, height = 1400, res = 300)
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
      just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  circlize_plot_MF(df_out_MF, df_genes_uni_MF, g)
  upViewport()
  draw(lgd_MF, x = circle_size, just = c("left", "top"))
  dev.off()
}


### oppA
# https://jokergoo.github.io/circlize_book/book/legends.html
library(gridBase)

g <- "oppA"
output_full_filepath <- "circular_oppA_BP_all.png"
create_plot_BP(output_full_filepath, df_outliers_BP_oppa_all, df_genes_unique_BP_oppa, lgd_list_vertical_BP_oppa, g)

output_full_filepath <- "circular_oppA_CC_all.png"
create_plot_CC(output_full_filepath, df_outliers_CC_oppa_all, df_genes_unique_CC_oppa, lgd_list_vertical_CC_oppa, g)

output_full_filepath <- "circular_oppA_MF_all.png"
create_plot_MF(output_full_filepath, df_outliers_MF_oppa_all, df_genes_unique_MF_oppa, lgd_list_vertical_MF_oppa, g)


### glgb
# https://jokergoo.github.io/circlize_book/book/legends.html
library(gridBase)

g <- "glgB"
output_full_filepath <- "circular_glgB_BP_all.png"
create_plot_BP(output_full_filepath, df_outliers_BP_glgb_all, df_genes_unique_BP_glgb, lgd_list_vertical_BP_glgb, g)

output_full_filepath <- "circular_glgB_CC_all.png"
create_plot_CC(output_full_filepath, df_outliers_CC_glgb_all, df_genes_unique_CC_glgb, lgd_list_vertical_CC_glgb, g)

output_full_filepath <- "circular_glgB_MF_all.png"
create_plot_MF(output_full_filepath, df_outliers_MF_glgb_all, df_genes_unique_MF_glgb, lgd_list_vertical_MF_glgb, g)


