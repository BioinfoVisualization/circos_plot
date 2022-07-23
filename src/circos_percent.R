#!/usr/bin/env Rscript
library(circlize)
library(RColorBrewer)
library(scales)
args = commandArgs(trailingOnly=TRUE)
path_to_file = args[1]
out_file = args[2]
circos.clear()
# Import adjacency matrix
mat = read.csv(path_to_file,header=TRUE,row.names=1)
mat <- as.matrix(mat)
# Set color pallete
n <- length(union(rownames(mat), colnames(mat)))
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
whites = c('#F5F5F5','#F7F7F7','#F6E8C3','#FDE0EF','#E7D4E8','#E6F5D0','#D9F0D3')
grid_col <- setdiff(col_vec,whites)[1:n]
# Circos plot
pdf(out_file, width = 20,height = 20)
chordDiagram(mat, grid.col = grid_col, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = mm_h(3)))



# Percentage
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  
  circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3) # dotted line
  by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.2, 0.5)
  for(p in seq(by, 1, by = by)) {
    circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.1, 
                paste0(p*100, "%"), cex = 1, adj = c(0.5, 0), niceFacing = TRUE)
  }
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(-0.35, 0.5),
              cex=1.5)}, bg.border = NA) # here set bg.border to NA is important

dev.off()
