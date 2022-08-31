library("data.table")
library("liftOver")

setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()
setwd("base data")
setwd("eRNA-regulation")
setwd("m3_regulation-TF")
heRaTF = readRDS("Pancreas-TF.rds")
setwd(path)
setwd("base data")
dataABC = fread("AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv")
pancGenes = subset(dataABC, dataABC$CellType == "body_of_pancreas-ENCODE")
setwd("eRNA-target")
setwd("m4_target")
heRaGene = readRDS("Pancreas-Target.rds")
setwd(path)
setwd("input")
geneList = fread("alle Gene.gene", header = F)
heRaGene = subset(heRaGene, heRaGene$Symbol %in% geneList$V1)
pancGenes = subset(pancGenes, pancGenes$TargetGene %in% geneList$V1)

rangesPancGenes = makeGRangesFromDataFrame(pancGenes)
overlaps = findOverlaps(rangesPancGenes)
test = as.matrix(overlaps)
View(test)
#no overlaps in ABC data, just discrete ranges mapped to multiple Genes (CTRB1/CTRB2 and CEL/CELP) (inversion mutation in CTRB1/2)

rangesHeRaGene = data.frame(heRaGene$echr, heRaGene$estart, heRaGene$eend)
rangesHeRaGene = makeGRangesFromDataFrame(rangesHeRaGene)
overlaps2 = findOverlaps(rangesHeRaGene)
test2 = as.matrix(overlaps2)
View(test2)

#heRa data show multiple overlaps within the same gene

overlaps3 = findOverlaps(rangesHeRaGene, rangesPancGenes)
test3 = as.matrix(overlaps3)
View(test3)
heRaGene[as.numeric(test3[1,1]),]
pancGenes[as.numeric(test3[1,2]),]
heRaGene[as.numeric(test3[2,1]),]
pancGenes[as.numeric(test3[2,2]),]
heRaGene[as.numeric(test3[2,1]),]
pancGenes[as.numeric(test3[2,2]),]

#one overlap between heRa and ABC in SPINK1 and two overlaps in CFTR intergenic enhancers
#TODO: overlap größe