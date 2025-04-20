library(ComplexHeatmap)
celltypes <- read.table('cell_type_use.txt')$V1
color <- read.table('ref_palette.new.txt', comment.char = '')
palatte <- color$V2
names(palatte) <- color$V1


celltype_pal <- c(
  "BMem"  = "#FF5733", 
  "BNav"  = "#33FF57", 
  "CD4EM"        = "#3357FF", 
  "CD4Naive"          = "#FF33FF",
  "CD4Treg"       = "#FF33A1", 
  "CD8GZMH"      = "#33FFA1", 
  "CD8GZMK"       = "#A1FF33", 
  "CD8Naive"    = "#FFA133",
  "MAIT"        = "#5733FF", 
  "MonocM"   = "#A133FF", 
  "NKBright"   = "#33A1FF", 
  "NKDim"    = "#FF5733",
  "NoClaM" = "#33FF57"
)  


label_dict<-list("BMem"="BMem",
                 "BNav"="BNav",
                 "CD4EM"="CD4EM",
                 "CD4Naive"="CD4Naive",
                 "CD4Treg"="CD4Treg",
                 "CD8GZMH"="CD8GZMH",
                 "CD8GZMK"="CD8GZMK",
                 "CD8Naive"="CD8Naive",
                 "MAIT"="MAIT",
                 "MonocM"="MonocM",
                 "NKBright"="NKBright",
                 "NKDim"="NKDim",
                 "NoClaM"="NoClaM")

dat <- read.delim('trans_sQTL.tsv')



for (i in seq_along(label_dict)) {
 dat$cell_type <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat$cell_type)
 celltypes <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],celltypes)
 names(palatte) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(palatte))
}


setOrder.manual <- c("BMem", "BNav", "CD4EM", "CD4Naive", "CD4Treg", "CD8GZMH", "CD8GZMK", "CD8Naive", "MAIT", "MonocM", "NKBright", "NKDim", "NoClaM")

lt = list()
for (ct in setOrder.manual) {
  tmp <- subset(dat, cell_type %in% ct)$gene_id
  lt <- append(lt, list(tmp))
}

names(lt) <- setOrder.manual

m1 = make_comb_mat(lt)
m <- m1[comb_size(m1) >= 3]

ss = set_size(m)
cs = comb_size(m)

pdf('Fig.5a.pdf', width = 6, height = 5)
UpSet(m, 
      set_order = match(setOrder.manual, set_name(m)), 
      right_annotation = upset_right_annotation(m, gp = gpar(col = NA, fill = palatte[names(ss)])),
      comb_order = order(comb_degree(m), -cs), 
      top_annotation = upset_top_annotation(m, gp = gpar(col = NA, fill = c(rep('black', length(grep(1, comb_degree(m), invert = T))), palatte[set_name(m)]))),
      comb_col = '#00004d'
)
dev.off()