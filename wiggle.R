#building the wiggle master file
#in a seperate script, because RAM usage is expected to be high
library("data.table")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("org.Hs.eg.db")

setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()

setwd("input")
genes = fread("alle Gene.gene", header = F)
dataABC = fread("AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv")
preFilteredGenes = subset(dataABC, dataABC$CellType == "body_of_pancreas-ENCODE")
#gene deren sequenzen zu untersuchen sind

setwd(gsub("/input","", getwd()))

setwd("output")

ensList = NULL
for (i in genes$V1){
  temp = mapIds(org.Hs.eg.db, keys = i,column = "ENTREZID", keytype = "SYMBOL")
  ensList = rbind(ensList, temp[[1]])
}
conversion = select(TxDb.Hsapiens.UCSC.hg38.knownGene, column = "TXSTART", keys = ensList, keytype = "GENEID")
uniqueConversion = unique(conversion$GENEID)
chromosomeMapping = select(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "CDSCHROM", keys = ensList, keytype = "GENEID")
chromosomeMapping = subset(chromosomeMapping, !is.na(chromosomeMapping$CDSCHROM))
chromosomeMapping = subset(chromosomeMapping, !grepl("alt", chromosomeMapping[,2]))
#script fuer ein sammelfile fuer alle seqs
chromosomes = paste("chr", c(1:22), sep="")
chromosomes = c(chromosomes, "chrX", "chrY")
for (i in chromosomes){
  tempWiggle = fread(i, header = T)
  for (i in subset(chromosomeMapping, chromosomeMapping$CDSCHROM == i)$GENEID){
    tss = subset(conversion, conversion$GENEID == i)
    tss = tss[1,2]
    start = tss - (kbRange * 1000)
    end = tss + (kbRange * 1000)
  }
}