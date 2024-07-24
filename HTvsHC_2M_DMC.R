if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")

setwd("path_to_working_directory") ### set up for DMCG analysis###

library(methylKit)

HC_Y1 <- "HC_Y1.sam"
HC_Y2 <- "HC_Y2.sam"
HC_Y3 <- "HC_Y3.sam"
HT_Y1 <- "HT_Y1.sam"
HT_Y2 <- "HT_Y2.sam"
HT_Y3 <- "HT_Y3.sam"

files_bam <- list(HC_Y1, HC_Y2, HC_Y3, HT_Y1, HT_Y2, HT_Y3)
my.methRaw=processBismarkAln(files_bam,
                             sample.id = list("HC_Y1", "HC_Y2", "HC_Y3", "HT_Y1", "HT_Y2","HT_Y3"),
                             assembly = "mm10",
                             read.context = "CpG",
                             mincov = 10,
                             minqual = 20,
                             treatment = c(0,0,0,1,1,1),
                             save.context = c("CpG"),
                             save.folder = getwd()
                             )

HC_Y1_CpG <- "HC_Y1_CpG.txt"
HC_Y2_CpG <- "HC_Y2_CpG.txt"
HC_Y3_CpG <- "HC_Y3_CpG.txt"
HT_Y1_CpG <- "HT_Y1_CpG.txt"
HT_Y2_CpG <- "HT_Y2_CpG.txt"
HT_Y3_CpG <- "HT_Y3_CpG.txt"


file_2.list <- list(HC_Y1_CpG, HC_Y2_CpG, HC_Y3_CpG, HT_Y1_CpG, HT_Y2_CpG, HT_Y3_CpG)
myobj = methRead(file_2.list,
                 sample.id = list("HC_Y1_CpG", "HC_Y2_CpG", "HC_Y3_CpG", "HT_Y1_CpG", "HT_Y2_CpG", "HT_Y3_CpG"),
                 assembly = "mm10",
                 pipeline = "bismark",
                 context = "CpG",
                 resolution = "base",
                 treatment = c(0,0,0,1,1,1),
                 mincov =10
                )

getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
getMethylationStats(myobj[[2]], plot = TRUE, both.strands = FALSE)
getMethylationStats(myobj[[3]], plot = TRUE, both.strands = FALSE)
getMethylationStats(myobj[[4]], plot = TRUE, both.strands = FALSE)
getMethylationStats(myobj[[5]], plot = TRUE, both.strands = FALSE)
getMethylationStats(myobj[[6]], plot = TRUE, both.strands = FALSE)

getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[2]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[3]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[4]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[5]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[6]], plot = TRUE, both.strands = FALSE)

filtered.myobj=filterByCoverage(myobj, 
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count = NULL,
                                hi.perc = 99.9
                               )

meth=unite(myobj, destrand = FALSE)

myDiff=calculateDiffMeth(meth)

myDiff25p=getMethylDiff(myDiff, difference = 25, qvalue = 0.05)
myDiff25p.hyper= getMethylDiff(myDiff, difference = 25, qvalue = 0.05, type = "hyper")
myDiff25p.hypo = getMethylDiff(myDiff, difference = 25, qvalue = 0.05, type = "hypo")
write.csv(myDiff, file = "myDiff.csv", row.names = FALSE)
write.csv(myDiff25p, file = "myDiff25p.csv", row.names = FALSE)
write.csv(myDiff25p.hyper, file = "myDiff25p_hyper.csv", row.names = FALSE)
write.csv(myDiff25p.hypo, file="myDiff25p_hypo.csv", row.names = FALSE)

#### methylKit offers methylation context string, ex: CpG,CHG,CHH, etc. (default:CpG)#####
#### therefore, set the path for DMCHG analysis as a new working directory for easiness #####

setwd("path_to_working_directory") 

my.methRaw=processBismarkAln(files_bam,
                             sample.id = list("HC_Y1", "HC_Y2", "HC_Y3", "HT_Y1", "HT_Y2","HT_Y3"),
                             assembly = "mm10",
                             read.context = "CHG",
                             mincov = 10,
                             minqual = 20,
                             treatment = c(0,0,0,1,1,1),
                             save.context = c("CHG"),
                             save.folder = getwd()
                             )

HC_Y1_CHG <- "HC_Y1_CHG.txt"
HC_Y2_CHG <- "HC_Y2_CHG.txt"
HC_Y3_CHG <- "HC_Y3_CHG.txt"
HT_Y1_CHG <- "HT_Y1_CHG.txt"
HT_Y2_CHG <- "HT_Y2_CHG.txt"
HT_Y3_CHG <- "HT_Y3_CHG.txt"


file_1.list <- list(HC_Y1_CHG, HC_Y2_CHG, HC_Y3_CHG, HT_Y1_CHG, HT_Y2_CHG, HT_Y3_CHG)

myobj=methRead(file_1.list,
               sample.id=list("HC_Y1_CHG", "HC_Y2_CHG", "HC_Y3_CHG", "HT_Y1_CHG", "HT_Y2_CHG", "HT_Y3_CHG"),
               assembly="mm10",
               pipeline="bismark",
               context="CHG",
               resolution="base",
               treatment = c(0,0,0,1,1,1),
               mincov=10
              )

filtered.myobj=filterByCoverage(myobj,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count = NULL,
                                hi.perc = 99.9
                               )

meth=unite(myobj, destrand = FALSE)

myDiff=calculateDiffMeth(meth)

myDiff25p=getMethylDiff(myDiff, difference = 25, qvalue = 0.05)
myDiff25p.hyper= getMethylDiff(myDiff, difference = 25, qvalue = 0.05, type = "hyper")
myDiff25p.hypo = getMethylDiff(myDiff, difference = 25, qvalue = 0.05, type = "hypo")
write.csv(myDiff, file = "myDiff.csv", row.names = FALSE)
write.csv(myDiff25p, file = "myDiff25p.csv", row.names = FALSE)
write.csv(myDiff25p.hyper, file = "myDiff25p_hyper.csv", row.names = FALSE)
write.csv(myDiff25p.hypo, file="myDiff25p_hypo.csv", row.names = FALSE)

#### Again set the path for DMCHH analysis as a new working directory for easiness #####

setwd("path_to_working_directory") 

my.methRaw=processBismarkAln(files_bam,
                             sample.id = list("HC_Y1", "HC_Y2", "HC_Y3", "HT_Y1", "HT_Y2","HT_Y3"),
                             assembly = "mm10",
                             read.context = "CHH",
                             mincov = 10,
                             minqual = 20,
                             treatment = c(0,0,0,1,1,1),
                             save.context = c("CHH"),
                             save.folder = getwd()
                             )

HC_Y1_CHH <- "HC_Y1_CHH.txt"
HC_Y2_CHH <- "HC_Y2_CHH.txt"
HC_Y3_CHH <- "HC_Y3_CHH.txt"
HT_Y1_CHH <- "HT_Y1_CHH.txt"
HT_Y2_CHH <- "HT_Y2_CHH.txt"
HT_Y3_CHH <- "HT_Y3_CHH.txt"


file_1.list <- list(HC_Y1_CHH, HC_Y2_CHH, HC_Y3_CHH, HT_Y1_CHH, HT_Y2_CHH, HT_Y3_CHH)

myobj=methRead(file_1.list,
               sample.id=list("HC_Y1_CHH", "HC_Y2_CHH", "HC_Y3_CHH", "HT_Y1_CHH", "HT_Y2_CHH", "HT_Y3_CHH"),
               assembly="mm10",
               pipeline="bismark",
               context="CHH",
               resolution="base",
               treatment = c(0,0,0,1,1,1),
               mincov=10
              )

filtered.myobj=filterByCoverage(myobj,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count = NULL,
                                hi.perc = 99.9
                               )

meth=unite(myobj, destrand = FALSE)

myDiff=calculateDiffMeth(meth)

myDiff25p=getMethylDiff(myDiff, difference = 25, qvalue = 0.05)
myDiff25p.hyper= getMethylDiff(myDiff, difference = 25, qvalue = 0.05, type = "hyper")
myDiff25p.hypo = getMethylDiff(myDiff, difference = 25, qvalue = 0.05, type = "hypo")
write.csv(myDiff, file = "myDiff.csv", row.names = FALSE)
write.csv(myDiff25p, file = "myDiff25p.csv", row.names = FALSE)
write.csv(myDiff25p.hyper, file = "myDiff25p_hyper.csv", row.names = FALSE)
write.csv(myDiff25p.hypo, file="myDiff25p_hypo.csv", row.names = FALSE)


