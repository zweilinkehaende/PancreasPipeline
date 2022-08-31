#coordinates first
#libraries
library("data.table")
library("liftOver")
library("parallel")

#define path
setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()

#base data
setwd("base data")

dataABC = fread("AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv")
preFilteredGenes = subset(dataABC, dataABC$CellType == "body_of_pancreas-ENCODE")
#Datenbankoutput der ABC Max Methode
TFBSs = fread("Human_TF_MotifList_v_1.01.csv")
TFBSs = subset(TFBSs, TFBSs$`Best Motif(s)? (Figure 2A)`== "TRUE")
getex = fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")
tfTPM = subset(getex, getex$Description %in% TFBSs$`HGNC symbol`)
tfTranscribedInP = subset(tfTPM, tfTPM$Pancreas > 1)

seqCoords = NULL

prep = function(relevantGenes, sqFilter, trimBool, trimLength, sqRange){
  #requires libraries "data.table" and "liftOver"
  #requires global variables "dataABC", "tfTranscribedInP" and "genome"
  
  setwd("input")
  genes = fread(relevantGenes, header = FALSE)
  #gene deren sequenzen zu untersuchen sind
  filteredGenes = subset(preFilteredGenes, preFilteredGenes$TargetGene %in% genes$V1)
  whitelist = NULL
  if (grepl("a", sqFilter)){
    whitelist = subset(filteredGenes, (filteredGenes$class == "intergenic" | filteredGenes$class == "genic"))
  }
  if (grepl("b", sqFilter)){
    transient = subset(filteredGenes, (filteredGenes$class == "intergenic" | filteredGenes$class == "genic"))
    for (i in unique(transient$TargetGene)){
      whitelist = rbind(whitelist, subset(transient, transient$"ABC.Score" == max(subset(transient, transient$TargetGene == i)$"ABC.Score")))
    }
  }
  if (grepl("p", sqFilter)){
    whitelist = rbind(whitelist, filteredGenes = subset(filteredGenes, filteredGenes$class == "promoter"))
  }
  if (grepl("t", sqFilter)){
    whitelist = rbind(whitelist, filteredGenes = subset(filteredGenes, filteredGenes$class == "promoter"))
  }
  
  setwd(gsub("/input","", getwd()))
  
  #creating target directory
  setwd("output")
  outputDirName = paste(gsub(".gene", "", relevantGenes), sqFilter)
  if (trimBool == 1){
    outputDirName = paste(outputDirName, trimLength)
  }
  if (grepl("t", sqFilter)){
    outputDirName = paste(outputDirName, sqRange)
  }
  dir.create(outputDirName)
  setwd(outputDirName)
  
  #liftover
  
  
  for (i in whitelist$start){
    whitelist$start = whitelist$start -1 #enumeration shift
    whitelist$end = whitelist$end -1
    whitelist$TargetGeneTSS = whitelist$TargetGeneTSS -1
  }
  if (grepl("t", sqFilter)){
    for (i in whitelist$TargetGeneTSS){
      whitelist$start = whitelist$TargetGeneTSS
      whitelist$end = whitelist$TargetGeneTSS +1
    }
  }
  newFormat = makeGRangesFromDataFrame(whitelist)
  liftedSeqs = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
  
  fileName = outputDirName
  fileName = paste(fileName, ".fa", sep="")
  #script fuer ein sammelfile fuer alle seqs
  x = 1
  for (i in whitelist$name){
    start = start(liftedSeqs[[x]])
    end = end(liftedSeqs[[x]])
    chr = subset(whitelist, whitelist$name == i)$chr
    chr = unique(chr)
    start = start +1 #enumeration shifted back to 1-indexing
    end = end +1
    sqName = subset(whitelist, whitelist$name == i)$TargetGene
    sqName = sqName[1]
    if (grepl("t", sqFilter)){
      sqName = paste(sqName, "_range_", sqRange, sep = "")
    }
    if (trimBool == 1){
      sqName = paste(sqName, "_trim_", trimLength, sep = "")
    }
    #start() and end() of liftedSeqs[[x]] sometimes return multiple results, first one is used in such cases
    center = floor(mean(c(start, end)))
    if (trimBool == 1){
      if ((end-center)>trimLength){
        start = center - trimLength
        end = center + trimLength
      }
    }
    if (grepl("t", sqFilter)){
      end = start + sqRange *1000
      start = start - sqRange *1000
    }
    # seqCoords <<- rbind(seqCoords, c(chr, outputDirName, sqName, start, end)) #including outputDirName causes the script to crash at "relevante Gene T 0.8", maybe caused by rbind using spaces as delimiters
    seqCoords <<- rbind(seqCoords, c(chr, paste(outputDirName), sqName, start, end))
    x = x +1
  }
  setwd(gsub(paste("/", outputDirName, sep=""),"", getwd()))
  setwd(gsub("/output","", getwd()))
  print(outputDirName)
}
setwd(path)
setwd("input")
inputList = list.files()
geneList = NULL
tfList = NULL
x = 1
y = 1
for (i in inputList){
  if (substring(i, nchar(i))=="e"){ #looking for .gene files
    geneList[x] = i
    x = x+1
  }
  if (substring(i, nchar(i))=="f"){ #looking for .tf files
    tfList[y] = i
    y = y+1
  }
}
setwd(gsub("/input","", getwd()))
regionsList = c(c(1:9)/10, c(1:9), c(1:10)*10)
trimList = c(250, 500, 1000)
sqFilterList = c("pa", "pb", "p", "a","b")
wiggleTable = NULL
seqCoords = NULL

setwd(path)
for (i in geneList){
  for (k in trimList){
    for (l in sqFilterList){
      prep(i, l, 1, k, 0)
    }
  }
  for (m in sqFilterList){
    prep(i, m, 0, 10000, 0)
  }
  
  for (j in regionsList){
    prep(i, "t", 0, 0, j)
  }
}
seqCoords = seqCoords[order(seqCoords[,1], seqCoords[,4]),]
extremesMatrix = NULL

for (i in unique(seqCoords[,1])){
  temp = subset(seqCoords, seqCoords[,1] == i)
  minimum = min(temp[,4])
  maximum = max(temp[,5])
  extremesMatrix <<- rbind(extremesMatrix, c(i, minimum, maximum))
  print(i)
}
print(extremesMatrix)
extremesMatrix = extremesMatrix[order(extremesMatrix[,1]),]
# wiggleFile = NULL
# setwd(path)
# setwd("base data")
# setwd("custom wiggles")
# files = list.files()
# for (i in files){
#   wiggleFile = rbind(wiggleFile, fread(i, header = F))
# }
# lengths = NULL
# setwd(path)
# setwd("base data")
# setwd("phyloP100way")
# setwd("chromosomes")
# wiggleFile = NULL
# for (i in unique(extremesMatrix[,1])){
#   chrWiggle = fread(i, header = T)
#   temp = subset(extremesMatrix, extremesMatrix[,1] == i)
#   tempSq = chrWiggle[temp[2]:temp[3]]
#   wiggleFile <<- rbind(wiggleFile, tempSq, use.names = F, col.names = F, row.names = F)
# }
# remove(chrWiggle)
# remove(dataABC)
# remove(getex)
# remove(getexMeans)
# remove(preFilteredGenes)
# gc()
# setwd(path)
# setwd("output")
# backup = wiggleFile
# test2 = as.data.frame(backup)
# wiggleFile = as.data.frame(wiggleFile)
# wiggleFile = as.data.frame(sapply(wiggleFile, as.numeric))
# test = subset(test2, is.na(wiggleFile))
# naCoords = which(is.na(wiggleFile))
# wiggleFile[wiggleFile == FALSE] = 0
# for (i in naCoords){ #way too slow, would take ~6 hours
#   wiggleFile[i,1] = paste(backup[i])
#   print(i)
# }
# x = 1
# while (x <= length(wiggleFile[,1])){
#   if (!(is.na(as.numeric(wiggleFile[x,1])))){
#     wiggleFile[x,1] = as.numeric(wiggleFile[x,1])
#   }
#   print(x)
#   x = x +1
# }
# write("fixedStep chrom=chr1 start=10918 step=1", "wiggleTest.wig", append = F)
# write(wiggleFile, "wiggleTesttxt.wig", append = F)
# write.table(wiggleFile, "wiggleTest.wig", col.names= F, row.names = F, append = F, quote = F)

x = 1
while (x <= length(extremesMatrix[,1])){
  lengths = rbind(lengths, as.numeric(extremesMatrix[x,3]) - as.numeric(extremesMatrix[x,2]))
  x = x +1
}
consequtiveLengths = 0
x = 1
while (x < length(lengths)){
  consequtiveLengths = rbind(consequtiveLengths, lengths[x])
  x = x + 1
}
x = 1
while (x <= length(extremesMatrix[,1])){
  extremesMatrix[x,2] = as.numeric(extremesMatrix[x,2]) - as.numeric(consequtiveLengths[x])
  extremesMatrix[x,3] = as.numeric(extremesMatrix[x,3]) - as.numeric(consequtiveLengths[x])
  x = x +1
}

# i = 1
# while (i <= length(test[,3])){
#   temp = subset(extremesMatrix, extremesMatrix[,1] == test[i,1])
#   temp = temp[,2]
#   test[i,4] = as.numeric(test[i,4]) - as.numeric(temp)
#   test[i,5] = as.numeric(test[i,5]) - as.numeric(temp)
#   i = i +1
# }
setwd(path)
setwd("base data")
setwd("hg38/FASTA")
chromosomes = paste("chr", c(1:22), sep="")
chromosomes = c(chromosomes, "chrX", "chrY")
genome = NULL
for (i in chromosomes){
  genome[i] = fread(paste(i,".fna", sep=""))
}

setwd(gsub("/hg38/FASTA","", getwd()))
setwd(gsub("/base data","", getwd()))
setwd("output")
x = 1
while (x <= length(seqCoords[,1])){
  wd = seqCoords[x,2]
  setwd(path)
  setwd("base data")
  setwd("phyloP100way")
  setwd("chromosomes")
  start = as.numeric(seqCoords[x,4])
  end = as.numeric(seqCoords[x,5])
  startRow = floor(start/80)+1 #+1 wegen floor()
  startPoint = start - startRow*80 + 80
  endRow = floor(end/80)+1
  endPoint = end - endRow*80 + 80
  fileName = paste(seqCoords[x, 2], ".fa", sep ="")
  targetRows = paste(unlist(eval(parse(text = paste("genome$", seqCoords[x,1], "[",startRow,":",endRow,"]", sep="")))), collapse="")
  targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
  setwd(path)
  setwd("output")
  setwd(wd)
  write(paste("./bigWigToWig -chrom=", seqCoords[x,1], " -start=", (start -1), " -end=", end, ' "base data/phyloP100way/hg38.phyloP100way.bw" "output/', wd, "/", start,'.wig"', sep = ""), "wigScript.sh", append = T)
  write(paste(">", seqCoords[x, 1], ":", seqCoords[x,4], "-", seqCoords[x,5], sep=""), file= fileName, append = T)
  write(targetSequence, file= fileName, append = T)
  write("", file=fileName, append = T)
  setwd("..")
  x = x +1
}
#Freeing up RAM (doesn't always work though)
rm(genome)
gc()
fileList = list.files()
setwd("..")

scriptStr = NULL
for (i in fileList){
  scriptStr = paste(scriptStr, './create-priors --oc "output/' , i, '" --parse-genomic-coord "', "output/", i, "/", i, '.fa" ', '"input/wiggleTest.wig"','\n', sep = "")
}
setwd("scripts")
threads = detectCores() -1
splitScripts = NULL
scriptStr = unlist(strsplit(scriptStr, "\n"))
x = 0
scriptParts = ceiling(length(scriptStr)/(threads)) 
write("#!/bin/bash\n", "createPriorsScriptOfScripts.sh")
while (x < (threads)){
  temp = eval(parse(text=paste("scriptStr[",x*scriptParts,":",(x+1)*scriptParts,"]", sep = "")))
  temp = temp[!is.na(temp)]
  write("#!/bin/bash", paste("createPriorsScriptPart", x, ".sh", sep = ""))
  write(temp, paste("createPriorsScriptPart", x, ".sh", sep = ""), append = T)
  write(paste("bash createPriorsScriptPart", x, ".sh", " & ", sep = ""), "createPriorsScriptOfScripts.sh", append = T)
  x = x + 1
}
setwd("..")
#TODO: Alles mit wiggle neu. wiggle files müssen mit der utility erstellt werden (script), damit diese jeweils die korrekten header haben
#TODO: wiggle files müssen für jedes FASTA-file seperat erstellt werden, weil create-priors streikt wenn mehr Daten als sequenzen angegeben werden.
#TODO: (einige einträge in seqCoords mit dem gleichen outputDirName sind duplikate (TSS)) Remove TSS multiplets
#TODO: Identifier in den sqName ob promotor oder enhancer
#DONE: Name der FASTA-seq in der dritten Spalte einfügen
#DONE: FASTA-files bauen
#DONE: wiggle-file bauen
#TODO: scripts bauen
