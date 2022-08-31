library("data.table")
library("liftOver")
library("parallel")

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
setwd("base data")
setwd("hg38/FASTA")
chromosomes = paste("chr", c(1:22), sep="")
chromosomes = c(chromosomes, "chrX", "chrY")
genome = NULL
for (i in chromosomes){
  genome[i] = fread(paste(i,".fna", sep=""))
}
setwd(path)
setwd("output_mcast")

mcast_ABC = function(geneList){
  # geneList = "A.gene"
  setwd(path)
  setwd("input")
  genes = fread(geneList, header = F)
  setwd(path)
  setwd("output_mcast")
  # gene = "CLDN2"
  dirName = gsub(".gene", "", geneList)
  dir.create(dirName)
  for (j in genes$V1){
    
    temp = NULL
    temp = subset(relevantGenes, relevantGenes$TargetGene == j)
    if (is.na(temp[1,1])){
      setwd(dirName)
      write(paste("no ABC data for ", j," with selected filters"), "warning.msg", append = T)
      setwd("..")
    }
    if (!(is.na(temp[1,1]))){
      # temp2 = subset(heRaTF, heRaTF$eRNA %in% temp$eRNA)
      temp$start = temp$start -1
      temp$end = temp$end -1
      x = 1
      ###ranges
      for (i in regionsList){
        temp$start[1] = temp$TargetGeneTSS[1] -1
        temp$end[1] = temp$TargetGeneTSS[1]
        newFormat = makeGRangesFromDataFrame(temp[1])
        newFormat = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
        setwd(dirName)
        start = min(start(newFormat)) +1
        end = max(end(newFormat)) +1
        end = start[[1]] + (i*1000)
        start = start[[1]] - (i*1000)
        identifier = paste("range_", i, sep="")
        fileName = paste(dirName, "_", identifier, ".fa", sep="")
        startRow = floor(start/80)+1 #+1 wegen floor()
        startPoint = start - startRow*80 + 80
        endRow = floor(end/80)+1
        endPoint = end - endRow*80 + 80
        chr = temp$echr
        write(paste(">", j, " ", temp$chr[1], ":", start, "-", end, sep=""), fileName, append=T)
        targetRows = paste(unlist(eval(parse(text = paste("genome$", temp$chr[1],"[",startRow,":",endRow,"]", sep="")))), collapse="")
        targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
        write(targetSequence, fileName, append=T)
        x = x +1
        setwd("..")
      }
    }
  }
}
setwd(path)
setwd("input")
geneListList = subset(list.files(), grepl(".gene", list.files()))
relevantGenes = pancGenes
regionsList = c(c(1:9)/10, c(1:9), c(1:10)*10)
for (i in geneListList){
  mcast_ABC(i)
}

setwd(path)
setwd("output_mcast")
dirList = list.files()
scriptStr = NULL

for (i in dirList){
  setwd(i)
  fileList = list.files()
  for (j in fileList){
    scriptStr = paste(scriptStr, 'mcast --synth --oc "output_mcast/', i, '/mcast_out_', j, '" "common files/besonders wichtige TFs.meme" "output_mcast/', i, "/", j,  '"\n', sep ="")
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
write("#!/bin/bash\n", "mcastScriptOfScripts.sh")
while (x < (threads)){
  temp = eval(parse(text=paste("scriptStr[",x*scriptParts,":",(x+1)*scriptParts,"]", sep = "")))
  temp = temp[!is.na(temp)]
  write("#!/bin/bash", paste("mcastScriptPart", x, ".sh", sep = ""))
  write(temp, paste("mcastScriptPart", x, ".sh", sep = ""), append = T)
  write(paste("bash mcastScriptPart", x, ".sh", " & ", sep = ""), "mcastScriptOfScripts.sh", append = T)
  x = x + 1
}
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