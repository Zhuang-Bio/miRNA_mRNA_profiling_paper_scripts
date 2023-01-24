#create link for the miRNAs
miRBase_LINK <- "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="
ensembl_LINK <- "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="

createLink <- function(val,link,label) {
  sprintf('<a href="%s%s" target="_blank">%s</a>',link,val,label)
}

#this function to draw heatmap using pheatmap
DEcomplexheatmap <- function(matZval, columnName, rowName, columnClu, rowClu, annoColumn) {
  out <-  pheatmap::pheatmap(
      matZval, #matrix data
      color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), #or you can set any color you like, color = colorRampPalette(c("navy","white","firebrick"))(90)
      #fontsize=13, 
      scale="none",  #character indicating if the values should be centered and scaled in either the row direction or the column direction, or none.
      show_colnames = columnName,
      show_rownames = rowName,  
      cluster_cols = columnClu, 
      #clustering_distance_cols = "euclidean", 
      cluster_rows = rowClu, 
      #clustering_distance_rows = "euclidean",
      clustering_method = "complete", 
      cellheight = 20, 
      annotation_col = annoColumn#, 
      #annotation_row = FALSE, 
      #annotation_colors = annoCol
  )
  return(out)
}

#heatmap using Heatmap from ComplexHeatmap package
DEmRNAheatmap_f <- function(matZval, columnName, rowName, columnClu, rowClu, annoColumn) {
  out <-  Heatmap(
      matZval, #Z-score matrix data
      color = colorRamp2(c(-2,0,2),c("#2171b5","#F7F7F7","#B2182B")), #or you can set any color you like, color = colorRampPalette(c("navy","white","firebrick"))(90)
      #fontsize=13, 
      name="Z-score",  #character indicating if the values should be centered and scaled in either the row direction or the column direction, or none.
      show_column_names = columnName,
      show_row_names = rowName,
      cluster_columns = columnClu,
      cluster_rows = rowClu,
      top_annotation = annoColumn,
      use_raster = TRUE
  )
  return(out)
}

vioplotFun <- function(datadf,aesX, aesY,compares){
  out <-  ggviolin(datadf, x = aesX, y = aesY, fill = aesX,
            palette = c('#e41a1c','#377eb8','#4daf4a','#984ea3'),
            add = "boxplot", add.params = list(fill = "white")) +
            stat_compare_means(comparisons = compares, label = "p.signif")
  return(out)
}

#htmapIndmiRsep <- function(trait, data, rowano, text1, text2) {
#  outht <- Heatmap(as.matrix(data), name = "-log10(Pvalue)",
#                          show_column_names = T, show_row_names = T,
#                          column_title = paste0(trait, ": Individual-miRNAs and mRNA-modules enrichments"),
#                          row_names_side = "right",
#                          column_names_rot = 45,
#                          left_annotation = rowano,
#                          cluster_columns = FALSE, 
#                          cluster_rows = FALSE,
#                          cell_fun = function(j, i, x, y, width, height, fill) {
#                            if(text1[i, j] > 1 & text2[i, j] > -log10(0.05))
#                              grid.text(sprintf("%s", text1[i, j]), x, y, gp = gpar(fontsize = 8))
#                          },
#                          col= colorRamp2(c(0,2,4,6),c("white","#ffd8d8","#ff7676","red")),
#                          border = TRUE)
#  return(outht)
#}



