#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(circlize))
library(circlize)
library(RColorBrewer)
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
whites = c('#F5F5F5','#F7F7F7','#F6E8C3','#FDE0EF','#E7D4E8','#E6F5D0')
grid_col <- setdiff(col_vec,whites)[1:n]
# Circos plot
pdf(out_file, width = 20,height = 20)
chordDiagram(mat,grid.col = grid_col, annotationTrack = "grid",
             annotationTrackHeight = mm_h(2),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
# Twist labels perpendicular to plot
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex=1.2)}, bg.border = NA) # here set bg.border to NA is important
dev.off()
print('Circos plot created')
