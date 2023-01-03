#!/usr/bin/env Rscript
f <- commandArgs(trailingOnly = TRUE)
# testing
#f <- "/mnt/ceph/454_data/dante/Intact_LTR_dataset/assembly_o_sativa/Oryza_sativa_msu7_dante_ltr_bugfix32_trna_hemi_27.gff3"

library(rtracklayer)
gff <- import(f, format = "gff3")
eval <- as.numeric(gff$evalue)

gff_te <- gff[gff$type == "transposable_element" & !grepl("partial", gff$ID) & !is.na(eval),]

is_hemi <- function(x) {
  # count lower case actg:
    lower_case <- sum(grepl("[actg]", strsplit(x, "")[[1]]))
    upper_case <- sum(grepl("[ACTG]", strsplit(x, "")[[1]]))
    if (lower_case / upper_case > 1){
        return(TRUE)
    } else {
        return(FALSE)
    }
}
# vectorize is_hemi
is_hemi <- Vectorize(is_hemi)


all_pbs <- as.matrix(table(gff_te$trna_id, gff_te$Final_Classification))
# sort by AA
AA <-  sapply(strsplit(sapply(strsplit(rownames(all_pbs),"__"),"[[", 2),"-"), "[[", 1)
AAA <-  sapply(strsplit(sapply(strsplit(gff_te$trna_id,"__"),"[[", 2),"-"), "[[", 1)
trna_type <- ifelse(is_hemi(rownames(all_pbs)), "hemi", "full")
all_pbs_sorted <- all_pbs[order(trna_type, AA),]
# export

write.table(all_pbs_sorted, file = paste0(f, "_pbs_all.csv"),
            sep = "\t", quote = FALSE, row.names = TRUE,
            col.names = TRUE)



gff_te_p <- split(gff_te, gff_te$Final_Classification)


pdf(paste0(f, "_pbs.pdf"), width = 11, height = 20)
for (n in names(gff_te_p)){
  pbs <- gff_te_p[[n]]$trna_id
  eval <- as.numeric(gff_te_p[[n]]$evalue)
  AAp <-  sapply(strsplit(sapply(strsplit(pbs,"__"),"[[", 2),"-"), "[[", 1)
  AA_levels <- names(sort(table(AAp)))
  if (all(eval >= 0.1)){
    AAp01 <- character(0)
  }else{
    AAp01 <-  sapply(strsplit(sapply(strsplit(pbs[eval < 0.1],"__"),"[[", 2),"-"), "[[", 1)
  }
  if (all(eval >= 0.5)){
    AAp05 <- character(0)
  }else{
    AAp05 <-  sapply(strsplit(sapply(strsplit(pbs[eval < 0.5],"__"),"[[", 2),"-"), "[[", 1)
  }
  par(mfrow = c(3,1), mar = c(4,20,2,2))
  barplot(sort(table(gff_te_p[[n]]$trna_id), decreasing = TRUE), main = n, horiz = TRUE, las = 1)
  AAp_all <- table(factor(AAp, levels = AA_levels))
  AAp_01 <- table(factor(AAp01, levels = AA_levels))
  AAp_05 <- table(factor(AAp05, levels = AA_levels))
  AAp2 <- rbind(AAp_all, AAp_01)
  barplot(AAp_all, main = n, horiz = TRUE, las = 1)
  barplot(AAp_05, main = n, horiz = TRUE, las = 1, add = TRUE, col = "orange")
  barplot(AAp_01, main = n, horiz = TRUE, las = 1, add = TRUE, col = "red")
  legend("bottomright", legend = c("all", "evalue < 0.5", "evalue < 0.1"), fill = c("gray", "orange", "red"))
  hist(-log10(as.numeric(gff_te_p[[n]]$evalue)), xlab = "-log10 evlaue", main = n, breaks = 100)
}
par(mfrow = c(3,1))
hist(-log10(as.numeric(gff_te$evalue)), xlab = "-log10 evlaue", main = "all", breaks = 20)
barplot(sort(table(AAA), decreasing = TRUE), main = "all", horiz = TRUE, las = 1)
dev.off()