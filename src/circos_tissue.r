
library(circlize)
library(RColorBrewer)

path_to_file = '../data/output_data/circos_adj_in_tissue/circos_adjacency_A_thresh_2.csv'
out_file = '../figures/circos_tissue/circos1.pdf'

circos.clear()
# Import adjacency matrix
mat = read.csv(path_to_file,header=TRUE,row.names=1)
mat <- as.matrix(mat)

# Set sparse interactions
c = ncol(mat)
r = nrow(mat)
link.visible = matrix(1, r, c)
rownames(link.visible) = rownames(mat)
colnames(link.visible) = colnames(mat)
link.visible[r,]=0
link.visible[,c]=0

# Set color pallete
n <- length(union(rownames(mat), colnames(mat)))
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
whites = c('#F5F5F5','#C7EAE5','#F7F7F7','#F6E8C3','#FDE0EF','#E7D4E8','#E6F5D0','#D9F0D3')
grid_col <- setdiff(col_vec,whites)[1:n]

# Circos plot
pdf(out_file, width = 20,height = 20)
circos.par(start.degree = -60, clock.wise = TRUE)
chordDiagram(mat,grid.col = grid_col, annotationTrack = "grid", big.gap = 2,
             link.visible = link.visible, annotationTrackHeight = mm_h(2),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
# Twist labels perpendicular to plot
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex=1.2)}, bg.border = NA) # here set bg.border to NA is important
dev.off()
print('Circos plot created')
