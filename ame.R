library("data.table")
library("liftOver")
library("parallel")

setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()
setwd(path)
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
setwd("base data")
setwd("hg38/FASTA")
chromosomes = paste("chr", c(1:22), sep="")
chromosomes = c(chromosomes, "chrX", "chrY")
genome = NULL
for (i in chromosomes){
  genome[i] = fread(paste(i,".fna", sep=""))
}
setwd(path)
setwd("output_ame")
ame_HeRA = function(gene){
  setwd(path)
  setwd("output_ame")
  # gene = "PRSS1"
  temp = NULL
  temp = subset(heRaGene, heRaGene$Symbol == gene)
  # temp2 = subset(heRaTF, heRaTF$eRNA %in% temp$eRNA)
  newFormat = NULL
  newFormat = data.frame(temp$echr, temp$estart, temp$eend)
  newFormat$temp.estart = newFormat$temp.estart -1
  newFormat$temp.eend = newFormat$temp.eend -1
  x = 1
  while (x <= length(newFormat$temp.estart)){
    newFormat = makeGRangesFromDataFrame(newFormat)
    newFormat = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
    dir.create(gene)
    setwd(gene)
    start = min(start(newFormat)) +1 #merging sequences split into multiples in hg19>hg38-liftover, since original data is from hg38 and continuous
    end = max(end(newFormat)) +1
    start = start[[1]]
    end = end[[1]]
    fileName = paste(gene, "_", x, ".fa", sep="")
    startRow = floor(start/80)+1 #+1 wegen floor()
    startPoint = start - startRow*80 + 80
    endRow = floor(end/80)+1
    endPoint = end - endRow*80 + 80
    chr = temp$echr
    write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), fileName, append=F)
    targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$echr,"[",startRow,":",endRow,"]", sep="")))), collapse="")
    targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
    write(targetSequence, fileName, append=T)
    x = x +1
    print(start)
    print(end)
  }
}
ame_ABC = function(gene){
  setwd(path)
  setwd("output_ame_ABC")
  # gene = "CLDN2"
  dir.create(gene)
  temp = NULL
  temp = subset(relevantGenes, relevantGenes$TargetGene == gene)
  if (is.na(temp[1,1])){
    setwd(gene)
    write("no ABC data with selected filters", "warning.msg")
    setwd("..")
  }
  if (!(is.na(temp[1,1]))){
    # temp2 = subset(heRaTF, heRaTF$eRNA %in% temp$eRNA)
    temp$start = temp$start -1
    temp$end = temp$end -1
    bestE = max(subset(enhancers, enhancers$TargetGene == gene)$"ABC.Score")
    bestP = max(subset(promoters, promoters$TargetGene == gene)$"ABC.Score")
    x = 1
    while (x <= length(temp$start)){
      newFormat = makeGRangesFromDataFrame(temp[x])
      newFormat = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
      
      setwd(gene)
      start = min(start(newFormat)) +1
      end = max(end(newFormat)) +1
      start = start[[1]]
      end = end[[1]]
      if (temp$class[x] == "intergenic" | temp$class[x] == "genic"){
        identifier = "e"
      }
      if (temp$class[x] == "intergenic" | temp$class[x] == "genic"){
        if (temp$"ABC.Score"[x] == bestE){
          identifier = "bestE"
        }
      }
      if (temp$class[x] == "promoter"){
        identifier = "p"
      }
      if (temp$class[x] == "promoter"){
        if (temp$"ABC.Score"[x] == bestP){
          identifier = "bestP"
        }
      }
      fileName = paste(gene, "_", x, "_", identifier, ".fa", sep="")
      startRow = floor(start/80)+1 #+1 wegen floor()
      startPoint = start - startRow*80 + 80
      endRow = floor(end/80)+1
      endPoint = end - endRow*80 + 80
      write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), fileName, append=F)
      targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
      targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
      write(targetSequence, fileName, append=T)
      if (!(identifier == "bestP")){
        write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), paste(gene, "_all_E.fa", sep = ""), append=T)
        targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
        targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
        write(targetSequence, paste(gene, "_all_E.fa", sep = ""), append=T)
      }
      write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), paste(gene, "_all.fa", sep = ""), append=T)
      targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
      targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
      write(targetSequence, paste(gene, "_all.fa", sep = ""), append=T)
      setwd("..")
      for (i in names(geneGroups)){
        if (gene %in% eval(parse(text = paste("geneGroups$'", i, "'", sep = "")))){
          dirName = gsub(".gene", "", i)
          setwd(dirName)
          if (!(identifier == "bestP")){
            write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), paste(i, "_all_E.fa", sep = ""), append=T)
            targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
            targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
            write(targetSequence, paste(i, "_all_E.fa", sep = ""), append=T)
          }
          if (identifier == "bestP"){
            write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), paste(i, "_all_P.fa", sep = ""), append=T)
            targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
            targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
            write(targetSequence, paste(i, "_all_P.fa", sep = ""), append=T)
          }
          write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), paste(i, "_all.fa", sep = ""), append=T)
          targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
          targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
          write(targetSequence, paste(i, "_all.fa", sep = ""), append=T)
          setwd("..")
        }
      }
      x = x +1
    }
    ###ranges
    for (i in regionsList){
      temp$start[1] = temp$TargetGeneTSS[1] -1
      temp$end[1] = temp$TargetGeneTSS[1]
      newFormat = makeGRangesFromDataFrame(temp[1])
      newFormat = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
      setwd(gene)
      start = min(start(newFormat)) +1
      end = max(end(newFormat)) +1
      end = start[[1]] + (i*1000)
      start = start[[1]] - (i*1000)
      identifier = paste("range_", i, sep="")
      fileName = paste(gene, "_", x, "_", identifier, ".fa", sep="")
      startRow = floor(start/80)+1 #+1 wegen floor()
      startPoint = start - startRow*80 + 80
      endRow = floor(end/80)+1
      endPoint = end - endRow*80 + 80
      chr = temp$echr
      write(paste(">", temp$chr[1], ":", start, "-", end, sep=""), fileName, append=F)
      targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
      targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
      write(targetSequence, fileName, append=T)
      x = x +1
      setwd("..")
    }
  }
}
setwd(path)
setwd("input")
geneList = fread("alle Gene in heRa.gene", header = F)

for (i in geneList$V1){
  ame_HeRA(i)
}
setwd(path)
setwd("input")
geneList = fread("alle Gene.gene", header = F)
tfListRelevant = fread("relevante TFs.tf", header = F)
relevantGenes = subset(pancGenes, pancGenes$TargetGene %in% geneList$V1)
enhancers = subset(relevantGenes, (relevantGenes$class == "intergenic" | relevantGenes$class == "genic"))
promoters = subset(relevantGenes, relevantGenes$class == "promoter")
regionsList = c(c(1:9)/10, c(1:9), c(1:10)*10)
setwd(path)
setwd("input")
fileList = subset(list.files(), grepl(".gene", list.files()))
geneGroups = NULL
for (i in fileList){
  geneGroups[i] = fread(i, header = F)
}
setwd(path)
setwd("output_ame_ABC")
for (i in fileList){
  dir.create(gsub(".gene", "", i))
}
# for (j in geneList$V1){
#   gene = j
#   # print(gene)
#   for (i in names(geneGroups)){
#     
#     if (gene %in% eval(parse(text = paste("geneGroups$'", i, "'", sep = "")))){
#       dirName = gsub(".gene", "", i)
#       print(dirName)
#     }
#   }
# }

for (i in geneList$V1){
  ame_ABC(i)
}


setwd(path)
setwd("output_ame")
dirList = list.files()
scriptStr = NULL
scoringMethod = "fisher"
# scoringMethod = "ranksum"
scoringArgument = paste(" --method ", scoringMethod, " ", sep = "")


for (i in dirList){
  setwd(i)
  fileList = list.files()
  for (j in fileList){
    scriptStr = paste(scriptStr, 'ame --control --shuffle--', scoringArgument, '--evalue-report-threshold 401 --oc "output_ame/', i, '/ame_out_', j, '" "output_ame/', i, "/", j,  '" "common files/allTFs.meme" \n', sep ="")
  }
  setwd("..")
}
setwd(path)
setwd("output_ame_ABC")
dirList = list.files()


for (i in dirList){
  setwd(i)
  fileList = list.files()
  for (j in fileList){
    scriptStr = paste(scriptStr, 'ame --control --shuffle--', scoringArgument, '--evalue-report-threshold 401 --oc "output_ame_ABC/', i, '/ame_out_', j, '" "output_ame_ABC/', i, "/", j,  '" "common files/allTFs.meme" \n', sep ="")
  }
  setwd("..")
}
setwd(path)
setwd("scripts")
threads = detectCores() -1
splitScripts = NULL
scriptStr = unlist(strsplit(scriptStr, "\n"))
x = 0
scriptParts = ceiling(length(scriptStr)/(threads)) 
write("#!/bin/bash\n", "ameScriptOfScripts.sh")
while (x < (threads)){
  temp = eval(parse(text=paste("scriptStr[",x*scriptParts,":",(x+1)*scriptParts,"]", sep = "")))
  temp = temp[!is.na(temp)]
  write("#!/bin/bash", paste("ameScriptPart", x, ".sh", sep = ""))
  write(temp, paste("ameScriptPart", x, ".sh", sep = ""), append = T)
  write(paste("bash ameScriptPart", x, ".sh", " & ", sep = ""), "ameScriptOfScripts.sh", append = T)
  x = x + 1
}

setwd("output_ame_ABC")
pathAnalysis = getwd()

setwd(pathAnalysis)
fileList = list.files()

fileList[1]
pHits = NULL
y = 0
for (i in fileList){
  setwd(pathAnalysis)
  setwd(i)
  fileListInternal = list.dirs()
  for (j in fileListInternal){
    if (j == "."){
      next
    }
    setwd(j)
    temp = fread("ame.tsv", header = T)
    temp = subset(temp, temp$`adj_p-value` < 0.05)
    y = y +1
    pHits = rbind(pHits, c(temp$FASTA_max [1], length(temp$rank)))
    setwd("..")
  }
}
pHits = as.data.frame(pHits)
colnames(pHits) = c("seqs", "hits")
plot(hits ~ seqs, pHits)
setwd(path)
write.csv(pHits, "pHits.csv")

setwd(pathAnalysis)
fileList = list.files()
fileList[1]
pSum = NULL
y = 0
for (i in fileList){
  setwd(pathAnalysis)
  setwd(i)
  fileListInternal = list.dirs()
  for (j in fileListInternal){
    if (j == "."){
      next
    }
    setwd(j)
    temp = fread("ame.tsv", header = T)
    tempSum = 0
    for (k in temp$`adj_p-value`){
      tempSum = tempSum + 1/as.numeric(k)
    }
    pSum = rbind(pSum, c(temp$FASTA_max [1], log(tempSum)))
    setwd("..")
  }
}
pSum = as.data.frame(pSum)
colnames(pSum) = c("seqs", "pSUM")
plot(pSUM ~ seqs, pSum)
setwd(path)
write.csv(pSum, "pSUM.csv")

setwd(pathAnalysis)
fileList = list.files()
fileList[1]
foundSeqs = NULL
y = 0
for (i in fileList){
  setwd(pathAnalysis)
  setwd(i)
  fileListInternal = list.dirs()
  for (j in fileListInternal){
    if (j == "."){
      next
    }
    setwd(j)
    temp = fread("ame.tsv", header = T)
    temp = subset(temp, temp$motif_ID %in% tfListRelevant$V1)
    temp = subset(temp, temp$`adj_p-value` < 0.05)
    y = y +1
    x = 1
    while (x < length(temp$motif_ID)){
      foundSeqs = rbind(foundSeqs, c(i, j, temp$motif_ID[x], temp$`adj_p-value`[x], temp$`E-value`[x]))
      x = x +1
    }
    setwd("..")
  }
}
foundSeqs = as.data.frame((foundSeqs))
groups = subset(foundSeqs, (grepl("all", foundSeqs[,2])))
groups = subset(foundSeqs, (grepl("all", foundSeqs[,2])|grepl("range", foundSeqs[,2])))
counts = NULL
for (i in unique(groups[,1])){
  temp = subset(groups, groups[,1] == i)
  # temp = subset(groups, groups[,1] == "A")
  tempE = subset(temp, grepl("E", temp[,2]))
  tempP = subset(temp, grepl("P", temp[,2]))
  tempA = subset(temp, !(grepl("P", temp[,2])))
  tempA = subset(tempA, !(grepl("E", tempA[,2])))
  tempA = subset(tempA, !(grepl("range", tempA[,2])))
  tempR = subset(temp, grepl("range", temp[,2]))
  tempCountsE = length(tempE[,1])
  tempCountsP = length(tempP[,1])
  tempCountsA = length(tempA[,1])
  tempCountsR = length(tempR[,1])
  counts = rbind(counts, c(i, "A", tempCountsA))
  counts = rbind(counts, c(i, "P", tempCountsP))
  counts = rbind(counts, c(i, "E", tempCountsE))
  counts = rbind(counts, c(i, "R", tempCountsR))
}
counts = as.data.frame(counts)
library(ggplot2)
toPlotA = subset(counts, counts[,2] == "A")
toPlotA = toPlotA[,3]
colnames(counts) = c("group", "level", "hit_count")
counts$level = as.factor(counts$level)
is.factor(counts$level)
counts$hit_count = as.numeric(counts$hit_count)
is.numeric(counts$hit_count)
boxplot(hit_count ~ level, data=counts)

median(subset(counts$hit_count, counts$level == "A"))
mean(subset(counts$hit_count, counts$level == "A"))
median(subset(counts$hit_count, counts$level == "P"))
mean(subset(counts$hit_count, counts$level == "P"))
median(subset(counts$hit_count, counts$level == "E"))
mean(subset(counts$hit_count, counts$level == "E"))
median(subset(counts$hit_count, counts$level == "R"))
mean(subset(counts$hit_count, counts$level == "R"))
boxplot()

# ame --control --shuffle-- --evalue-report-threshold 100 "A 0.1.fa" "besonders wichtige TFs.meme"
# runtime including ranges and allTFs on 15 threads on 8 zen+ cores: ~9h all threads, >10h for some threads

#DONE: split the output to create a single fasta file for every sequence
#DONE: merge sequences split by liftover
#DONE: Annotate best predicted promotor
#DONE: check for overlaps between enhancers
#DONE: Do the same for the ABC data + additional identifiers (p/a/b/t/gene)
#DONE: build the AME script
#DONE: try different scoring methods, only fisher and ranksum applicable to single seqs without FASTA scores 
#DONE: ranges um TSS
#DONE: auf "allTFs.meme" umstellen