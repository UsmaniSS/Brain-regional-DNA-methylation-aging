if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")

setwd("path_to_working_directory") ### set up for DMCG analysis###

library(methylKit)

HC_O1 <- "HC_O1.sam"
HC_O2 <- "HC_O2.sam"
HC_O3 <- "HC_O3.sam"
HT_O1 <- "HT_O1.sam"
HT_O2 <- "HT_O2.sam"
HT_O3 <- "HT_O3.sam"

files_bam <- list(HC_O1, HC_O2, HC_O3, HT_O1, HT_O2, HT_O3)
my.methRaw=processBismarkAln(files_bam,
                             sample.id = list("HC_O1", "HC_O2", "HC_O3", "HT_O1", "HT_O2","HT_O3"),
                             assembly = "mm10",
                             read.context = "CpG",
                             mincov = 10,
                             minqual = 20,
                             treatment = c(0,0,0,1,1,1),
                             save.context = c("CpG"),
                             save.folder = getwd()
                             )

HC_O1_CpG <- "HC_O1_CpG.txt"
HC_O2_CpG <- "HC_O2_CpG.txt"
HC_O3_CpG <- "HC_O3_CpG.txt"
HT_O1_CpG <- "HT_O1_CpG.txt"
HT_O2_CpG <- "HT_O2_CpG.txt"
HT_O3_CpG <- "HT_O3_CpG.txt"


file_2.list <- list(HC_O1_CpG, HC_O2_CpG, HC_O3_CpG, HT_O1_CpG, HT_O2_CpG, HT_O3_CpG)
myobj = methRead(file_2.list,
                 sample.id = list("HC_O1_CpG", "HC_O2_CpG", "HC_O3_CpG", "HT_O1_CpG", "HT_O2_CpG", "HT_O3_CpG"),
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
                             sample.id = list("HC_O1", "HC_O2", "HC_O3", "HT_O1", "HT_O2","HT_O3"),
                             assembly = "mm10",
                             read.context = "CHG",
                             mincov = 10,
                             minqual = 20,
                             treatment = c(0,0,0,1,1,1),
                             save.context = c("CHG"),
                             save.folder = getwd()
                             )

HC_O1_CHG <- "HC_O1_CHG.txt"
HC_O2_CHG <- "HC_O2_CHG.txt"
HC_O3_CHG <- "HC_O3_CHG.txt"
HT_O1_CHG <- "HT_O1_CHG.txt"
HT_O2_CHG <- "HT_O2_CHG.txt"
HT_O3_CHG <- "HT_O3_CHG.txt"


file_1.list <- list(HC_O1_CHG, HC_O2_CHG, HC_O3_CHG, HT_O1_CHG, HT_O2_CHG, HT_O3_CHG)

myobj=methRead(file_1.list,
               sample.id=list("HC_O1_CHG", "HC_O2_CHG", "HC_O3_CHG", "HT_O1_CHG", "HT_O2_CHG", "HT_O3_CHG"),
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
                             sample.id = list("HC_O1", "HC_O2", "HC_O3", "HT_O1", "HT_O2","HT_O3"),
                             assembly = "mm10",
                             read.context = "CHH",
                             mincov = 10,
                             minqual = 20,
                             treatment = c(0,0,0,1,1,1),
                             save.context = c("CHH"),
                             save.folder = getwd()
                             )

HC_O1_CHH <- "HC_O1_CHH.txt"
HC_O2_CHH <- "HC_O2_CHH.txt"
HC_O3_CHH <- "HC_O3_CHH.txt"
HT_O1_CHH <- "HT_O1_CHH.txt"
HT_O2_CHH <- "HT_O2_CHH.txt"
HT_O3_CHH <- "HT_O3_CHH.txt"


file_1.list <- list(HC_O1_CHH, HC_O2_CHH, HC_O3_CHH, HT_O1_CHH, HT_O2_CHH, HT_O3_CHH)

myobj=methRead(file_1.list,
               sample.id=list("HC_O1_CHH", "HC_O2_CHH", "HC_O3_CHH", "HT_O1_CHH", "HT_O2_CHH", "HT_O3_CHH"),
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
