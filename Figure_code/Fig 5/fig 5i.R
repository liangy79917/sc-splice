
pdf(paste(outdir, 'circos_plot.pdf', sep=""), width = 8, height = 8)
cisgene <- cisgene %>% select(chr, start, end, label, ppo)
transgene <- transgene %>% select(chr, start, end, label, ppo)
print(cisgene)
print(transgene)
circos.initializeWithIdeogram(plotType = "ideogram")

chromosomes <- paste0("chr", 1:22)
circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA, 
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    if (chr %in% chromosomes) {
      xcenter = CELL_META$xcenter
      circos.text(
        x = xcenter, 
        y = 1, 
        labels = gsub("chr", "", chr),
        facing = "bending.outside", 
        cex = 0.7
      )
    }
  }
)

circos.genomicLink(cisgene, transgene, col = cell$color, ppo1 = cisgene$ppo, ppo2 = transgene$ppo)
circos.clear()
dev.off()
