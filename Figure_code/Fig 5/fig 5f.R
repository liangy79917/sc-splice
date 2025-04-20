
library(dplyr)
library(ggplot2)
library(data.table)
library(readr)
library(locuscomparer)
library(cowplot)

sub <- function(row) {
        chr <- as.character(row["VID"])
        chr <- substr(chr, 1, nchar(chr) - 2)
        return(chr)
}

add_label <- function(merged, snp){
    merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, '')
    return(merged)
}

total <- fread("rsid_ref.txt")
 
face <- data.frame(rsid = character(), pos = integer(), pval1 = double(), pval2 = double(), logp1 = double(), logp2 = double())

merged <- fread("merged_IVNS1ABP_RPS24.txt",sep="\t")
merged <- merge(total,merged,by.x = "Variation ID", by.y = "VID", all = FALSE)


ld <- fread("ld.tsv",sep="\t")

colnames(merged)[which(colnames(merged) == "dbSNP")] <- "rsid"
merged <- merged %>% select(rsid, pos, epval, spval, logpe, logps)
colnames(merged) <- c("rsid", "pos", "pval1", "pval2", "logp1", "logp2")
merged$chr <- "1"
merged <- merged[(!duplicated(merged[, c("pos")])), ]
legend <- TRUE
legend_position <- c('bottomright','topright','topleft')
snp <- merged[which.min(merged$pval1*merged$pval2), 'rsid']

snp = get_lead_snp(merged, snp)
color <- assign_color(merged$rsid, snp, ld)

shape <- ifelse(merged$rsid == snp, 23, 21)
names(shape) <- merged$rsid

size <- ifelse(merged$rsid == snp, 4, 3)
names(size) <- merged$rsid

merged <- add_label(merged, snp)

p <- make_scatterplot(merged, 'IVNS1ABP cis-eQTL', 'RPS24 trans-sQTL', color,
                      shape, size, legend, legend_position)
ggsave("single.pdf",p, device = "pdf", width = 4, height = 4)


p2 <- make_combined_plot(merged, 'IVNS1ABP cis-eQTL', 'RPS24 trans-sQTL', ld, chr=1, snp = "rs10798014", combine = FALSE, legend)

library(ggplot2)
ggsave("locuscompare.pdf", plot = p2$locuscompare, device = "pdf", width = 6, height = 6)
ggsave("locuszoom1.pdf", plot = p2$locuszoom1, device = "pdf", width = 6, height = 6)
ggsave("locuszoom2.pdf", plot = p2$locuszoom2, device = "pdf", width = 6, height = 6)
