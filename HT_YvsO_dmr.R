setwd("path_to_files")
library(methylKit)
HT_Y1_CpG <- "HT_Y1_CpG.txt"
HT_Y2_CpG <- "HT_Y2_CpG.txt"
HT_Y3_CpG <- "HT_Y3_CpG.txt"
HT_O1_CpG <- "HT_O1_CpG.txt"
HT_O2_CpG <- "HT_O2_CpG.txt"
HT_O3_CpG <- "HT_O3_CpG.txt"

file_1.list <- list(HT_Y1_CpG, HT_Y2_CpG, HT_Y3_CpG, HT_O1_CpG, HT_O2_CpG, HT_O3_CpG)
myobj=methRead(file_1.list, 
               sample.id=list("HT_Y1_CpG", "HT_Y2_CpG", "HT_Y3_CpG", "HT_O1_CpG", "HT_O2_CpG", "HT_O3_CpG"),
               assembly="mm10",
               pipeline="bismark", 
               context="CpG", 
               resolution="base",
               treatment = c(0,0,0,1,1,1),
               mincov=10)

filtered.myobj=filterByCoverage(myobj, 
                                lo.count=10, 
                                lo.perc=NULL,
                                hi.count = NULL, 
                                hi.perc = 99.9)

meth=unite(myobj, 
           destrand = FALSE)

myDiff=calculateDiffMeth(meth)
myDiff10p=getMethylDiff(myDiff, difference = 10, qvalue = 0.05)
myDiff10p.hyper= getMethylDiff(myDiff, difference = 10, qvalue = 0.05, type = "hyper")
myDiff10p.hypo = getMethylDiff(myDiff, difference = 10, qvalue = 0.05, type = "hypo")
write.csv(myDiff, file = "myDiff.csv", row.names = FALSE)
write.csv(myDiff10p, file = "myDiff10p.csv", row.names = FALSE)
write.csv(myDiff10p.hyper, file = "myDiff10p_hyper.csv", row.names = FALSE)
write.csv(myDiff10p.hypo, file="myDiff10p_hypo.csv", row.names = FALSE)

chr_wise_10p <- diffMethPerChr(myDiff, plot = FALSE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
write.csv(chr_wise_10p, file = "chr_wise_10p.csv", row.names = FALSE)

library(dplyr)

myobj %>%
  filterByCoverage(lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9) %>%
  tileMethylCounts(win.size = 1000, step.size = 1000, cov.bases = 0) ->
  tiles_CpG
meth_tiles_CpG = unite(tiles_CpG, destrand = FALSE)
myDiff_tiles_CpG=calculateDiffMeth(meth_tiles_CpG)
myDiff_tiles_CpG_10p=getMethylDiff(myDiff_tiles_CpG, difference = 10, qvalue = 0.05)
myDiff_tiles_CpG_10p.hyper=getMethylDiff(myDiff_tiles_CpG, difference = 10, qvalue = 0.05, type = "hyper")
myDiff_tiles_CpG_10p.hypo=getMethylDiff(myDiff_tiles_CpG, difference = 10, qvalue = 0.05, type = "hypo")
write.csv(myDiff_tiles_CpG_10p, file = "myDiff_tiles_CpG_10p.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_10p.hyper, file = "myDiff_tiles_CpG_10phyper.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_10p.hypo, file = "myDiff_tiles_CpG_10phypo.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG, file="myDiff_tiles_CpG.csv", row.names = FALSE)
diffMethPerChr(myDiff_tiles_CpG, plot = TRUE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
chr_wise_w1000_s1000_10p <- diffMethPerChr(myDiff_tiles_CpG, plot = FALSE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
write.csv(chr_wise_w1000_s1000_10p, file = "chr_wise_w1000_s1000_10p.csv", row.names = FALSE)

myobj %>%
  filterByCoverage(lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9) %>%
  tileMethylCounts(win.size = 500, step.size = 500, cov.bases = 0) ->
  tiles_CpG_w500
meth_tiles_CpG_w500 = unite(tiles_CpG_w500, destrand = FALSE)
myDiff_tiles_CpG_w500=calculateDiffMeth(meth_tiles_CpG_w500)
myDiff_tiles_CpG_w500_10p=getMethylDiff(myDiff_tiles_CpG_w500, difference = 10, qvalue = 0.05)
myDiff_tiles_CpG_w500_10p.hyper=getMethylDiff(myDiff_tiles_CpG_w500, difference = 10, qvalue = 0.05, type = "hyper")
myDiff_tiles_CpG_w500_10p.hypo=getMethylDiff(myDiff_tiles_CpG_w500, difference = 10, qvalue = 0.05, type = "hypo")
write.csv(myDiff_tiles_CpG_w500_10p, file = "myDiff_tiles_CpG_w500_10p.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_w500_10p.hyper, file = "myDiff_tiles_CpG_w500_10phyper.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_w500_10p.hypo, file = "myDiff_tiles_CpG_w500_10phypo.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_w500, file="myDiff_tiles_CpG_w500.csv", row.names = FALSE)
diffMethPerChr(myDiff_tiles_CpG_w500, plot = TRUE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
chr_wise_w500_10p <- diffMethPerChr(myDiff_tiles_CpG_w500, plot = FALSE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
write.csv(chr_wise_w500_10p, file = "chr_wise_w500_10p.csv", row.names = FALSE)


myobj %>%
  filterByCoverage(lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9) %>%
  tileMethylCounts(win.size = 300, step.size = 300, cov.bases = 0) ->
  tiles_CpG_w300
meth_tiles_CpG_w300 = unite(tiles_CpG_w300, destrand = FALSE)
myDiff_tiles_CpG_w300=calculateDiffMeth(meth_tiles_CpG_w300)
myDiff_tiles_CpG_w300_10p=getMethylDiff(myDiff_tiles_CpG_w300, difference = 10, qvalue = 0.05)
myDiff_tiles_CpG_w300_10p.hyper=getMethylDiff(myDiff_tiles_CpG_w300, difference = 10, qvalue = 0.05, type = "hyper")
myDiff_tiles_CpG_w300_10p.hypo=getMethylDiff(myDiff_tiles_CpG_w300, difference = 10, qvalue = 0.05, type = "hypo")
write.csv(myDiff_tiles_CpG_w300_10p, file = "myDiff_tiles_CpG_w300_10p.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_w300_10p.hyper, file = "myDiff_tiles_CpG_w300_10phyper.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_w300_10p.hypo, file = "myDiff_tiles_CpG_w300_10phypo.csv", row.names = FALSE)
write.csv(myDiff_tiles_CpG_w300, file="myDiff_tiles_CpG_w300.csv", row.names = FALSE)
diffMethPerChr(myDiff_tiles_CpG_w300, plot = TRUE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
chr_wise_w300_10p <- diffMethPerChr(myDiff_tiles_CpG_w300, plot = FALSE, qvalue.cutoff = 0.05, meth.cutoff = 10, keep.empty.chrom = TRUE)
write.csv(chr_wise_w300_10p, file = "chr_wise_w300_10p.csv", row.names = FALSE)
