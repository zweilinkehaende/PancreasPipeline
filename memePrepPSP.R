library("data.table")
library("liftOver")
# library("biomaRt")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("org.Hs.eg.db")
library(parallel)

# test <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# getGene("HNF1A", type = 'hgnc_symbol', mart)
setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()
setwd("base data")
# HGNCtoENS = fread("ens_hgnc")

dataABC = fread("AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv")
preFilteredGenes = subset(dataABC, dataABC$CellType == "body_of_pancreas-ENCODE")
#Datenbankoutput der ABC Max Methode
TFBSs = fread("Human_TF_MotifList_v_1.01.csv")
TFBSs = subset(TFBSs, TFBSs$`Best Motif(s)? (Figure 2A)`== "TRUE")
getex = fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")
tfTPM = subset(getex, getex$Description %in% TFBSs$`HGNC symbol`)
tfTranscribedInP = subset(tfTPM, tfTPM$Pancreas > 1)
setwd("hg38/FASTA")
chromosomes = paste("chr", c(1:22), sep="")
chromosomes = c(chromosomes, "chrX", "chrY")
genome = NULL
for (i in chromosomes){
  genome[i] = fread(paste(i,".fna", sep=""))
}

setwd(gsub("/hg38/FASTA","", getwd()))
setwd(gsub("/base data","", getwd()))
hg38Annotated = (TxDb.Hsapiens.UCSC.hg38.knownGene)

#defines function to create seq file

prepMcast = function(geneGroup, kbRange){
  #requires libraries "data.table" and "liftOver"
  #requires global variables "dataABC", "tfTranscribedInP" and "genome"
  
  setwd("input")
  geneListOfGroup = fread(geneGroup, header = FALSE)
  geneListOfGroup = fread("alle Gene.gene", header = FALSE)
  #gene deren sequenzen zu untersuchen sind
  
  setwd(gsub("/input","", getwd()))
  
  outputDirName = paste(gsub(".gene", "", geneGroup), kbRange)
  setwd("output")
  dir.create(outputDirName)
  setwd(outputDirName)
  fileName = paste(outputDirName, ".fa", sep = "")
  ensList = NULL
  for (i in geneListOfGroup$V1){
    temp = mapIds(org.Hs.eg.db, keys = i,column = "ENTREZID", keytype = "SYMBOL")
    ensList = rbind(ensList, temp[[1]])
  }
  conversion = select(TxDb.Hsapiens.UCSC.hg38.knownGene, column = "TXSTART", keys = ensList, keytype = "GENEID")
  uniqueConversion = unique(conversion$GENEID)
  chromosomeMapping = select(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "CDSCHROM", keys = ensList, keytype = "GENEID")
  chromosomeMapping = subset(chromosomeMapping, !is.na(chromosomeMapping$CDSCHROM))
  chromosomeMapping = subset(chromosomeMapping, !grepl("alt", chromosomeMapping[,2]))
  #script fuer ein sammelfile fuer alle seqs
  x = 1
  write("", fileName, append = F)
  for (i in uniqueConversion){
    tempChromosomeMapping = subset(chromosomeMapping, chromosomeMapping$GENEID == i)
    sqName = geneListOfGroup$V1[x]
    if (is.na(tempChromosomeMapping[1,2])){
      write(paste("Warning Gene ID", i, "maps to no chr"), "warning.txt")
      next
    }
    tss = subset(conversion, conversion$GENEID == i)
    tss = tss[1,2]
    start = tss - (kbRange * 1000)
    end = tss + (kbRange * 1000)
    startRow = floor(start/80)+1 #+1 wegen floor()
    startPoint = start - startRow*80 + 80
    endRow = floor(end/80)+1
    endPoint = end - endRow*80 + 80
    targetRows = paste(unlist(eval(parse(text = paste("genome$", tempChromosomeMapping[1, 2],"[",startRow,":",endRow,"]", sep="")))), collapse="")
    targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
    write(paste(">", sqName, " ", start," Homo sapiens ", tempChromosomeMapping[1, 2],", GRCh38.p14", sep=""), file= fileName, append = T)
    write(targetSequence, file= fileName, append = T)
    write("", file=fileName, append = T)
    wiggleTable <<- rbind(wiggleTable, c(sqName, tempChromosomeMapping[1, 2], kbRange, start, end))
    x = x+1
    
  }
  setwd("..")
  setwd(gsub("/output","", getwd()))
  print(outputDirName)
}
prep = function(relevantGenes, sqFilter, trimBool, trimLength){
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

  setwd(gsub("/input","", getwd()))
  
  #creating target directory
  setwd("output")
  outputDirName = paste(gsub(".gene", "", relevantGenes), sqFilter)
  if (trimBool == 1){
    outputDirName = paste(outputDirName, trimLength)
  }
  dir.create(outputDirName)
  setwd(outputDirName)
  
  #liftover
  

  for (i in whitelist$start){
    whitelist$start = whitelist$start -1 #enumeration shift
    whitelist$end = whitelist$end -1
  }
  shiftedGenes = whitelist
  
  newFormat = makeGRangesFromDataFrame(shiftedGenes)
  liftedSeqs = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
  
  fileName = outputDirName
  fileName = paste(fileName, ".fa", sep="")
  #script fuer ein sammelfile fuer alle seqs
  x = 1
  for (i in whitelist$name){
    start = start(liftedSeqs[[x]])
    end = end(liftedSeqs[[x]])
    sqName = subset(whitelist, whitelist$name == i)
    #start() and end() of liftedSeqs[[x]] sometimes return multiple results, first one is used in such cases
    center = floor(mean(c(start, end)))
    if (trimBool == 1){
      if ((end-center)>trimLength){
        start = center - trimLength
        end = center + trimLength
      }
    }
    
    startRow = floor(start/80)+1 #+1 wegen floor()
    startPoint = start - startRow*80 +80
    endRow = floor(end/80)+1
    endPoint = end - endRow*80 +80
    targetRows = paste(unlist(eval(parse(text = paste("genome$",subset(whitelist, whitelist$name ==i)$chr,"[",startRow,":",endRow,"]", sep="")))), collapse="")
    targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
    write(paste(">", sqName$TargetGene, " ", start," Homo sapiens ", subset(whitelist, whitelist$name ==i)$chr,", GRCh38.p14", sep=""), file= fileName, append = T)
    write(targetSequence, file= fileName, append = T)
    write("", file=fileName, append = T)
    temp = subset(whitelist, whitelist$name ==i)
    wiggleTable <<- rbind(wiggleTable, c(sqName$TargetGene, temp$chr, sqFilter, start, end))
    x = x+1
    
  }
  setwd(gsub(paste("/", outputDirName, sep=""),"", getwd()))
  setwd(gsub("/output","", getwd()))
  print(outputDirName)
}
############
#builds the shared .meme motif files
setwd("input")
tfList = list.files()
setwd("..")
tfFilterListe = subset(tfList, grepl(".tf", tfList))
for (i in tfFilterListe){
  setwd("input")
  temp = fread(i, header = F)
  setwd("..")
  if (!(i == "allTFs.txt")){
    hgncToCisbp = subset(TFBSs, TFBSs$`HGNC symbol` %in% temp$V1)
  }
  if (i == "allTFs.txt"){
    hgncToCisbp = TFBSs
  }
  hgncToCisbpList = subset(hgncToCisbp, hgncToCisbp$`HGNC symbol` %in% tfTranscribedInP$Description)$"CIS-BP ID"
  setwd("base data/PSPMs")
  fileList = list.files()
  fileList = subset(fileList, gsub(".txt","",fileList) %in% hgncToCisbpList)
  setwd("..")
  setwd("..")
  
  setwd("output")
  fileName = gsub(".tf", "", i)
  write("MEME version 4", file= paste(fileName,".meme", sep=""))
  write("ALPHABET= ACGT", file= paste(fileName,".meme", sep=""), append =T)
  setwd("..")
  
  for (j in fileList){
    #TODO: pruefen ob das jeweilige motif die bedingungen von SEA an eine PSPM erfuellt (keine probability = 1, keine summe != 1)
    setwd("base data/PSPMs")
    motifName = j
    motifName = gsub(".txt","", motifName)
    motif = fread(j)
    pwm = motif[,2:5]
    pwmNameless = unname(pwm)
    identifier = subset(hgncToCisbp, hgncToCisbp$`CIS-BP ID` == gsub(".txt","",j))$"HGNC symbol"[1]
    alternateName = motifName
    setwd(gsub("base data/PSPMs","",getwd()))
    
    setwd("output")
    write(paste("MOTIF", identifier, alternateName, sep=" "), file=paste(fileName, ".meme", sep=""), append=T)
    write("letter-probability matrix:", file=paste(fileName, ".meme", sep=""), append=T)
    write.table(pwm,file=paste(fileName, ".meme", sep=""), append= T, col.names= F, row.names = F)
    write("\n", file=paste(fileName, ".meme", sep=""), append=T)
    setwd("..")
  }
}
##############################
#build the directory structure and seq.gene files
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

setwd(path)
for (i in geneList){
  for (k in trimList){
    for (l in sqFilterList){
      prep(i, l, 1, k)
    }
  }
  for (m in sqFilterList){
    prep(i, m, 0, 10000)
  }

  for (j in regionsList){
    prepMcast(i, (j))
  }
}
write.table(wiggleTable, "wiggleTable.txt", col.names= F, row.names = F)

##writes the script to run sea for all combinations
setwd(path)
setwd("output")
getwd()
threads = detectCores() -1 #using 1 thread fewer than available to prevent the computer from freezing up
#careful about RAM usage, each thread used close to 2 GB in testing, highly dependent on motif count. less than 5 GB total with < 100 motifs
#runtime on stock AMD 2700X: ~7 mins all threads, ~4 mins until completion.
dirs = subset(list.files(), !(grepl(".meme", list.files())))
scriptStr = NULL
for (i in dirs){
  sharedFiles = subset(list.files(), grepl(".meme", list.files()))
  
  setwd(paste(i))
  files = list.files()
  scriptStr = paste(scriptStr, 'mcast ','--synth --oc "mcast out/' , i, '" "besonders wichtige TFs.meme" "' , i, "/", i, '.fa"', '\n', sep = "")
  setwd("..")
  getwd()
}
setwd("..")
getwd()
setwd("scripts")
getwd()
splitScripts = NULL
scriptStr = unlist(strsplit(scriptStr, "\n"))
x = 0
scriptParts = ceiling(length(scriptStr)/(threads)) 
write("#!/bin/bash\n", "scriptOfScripts.sh")
while (x < (threads)){
  temp = eval(parse(text=paste("scriptStr[",x*scriptParts,":",(x+1)*scriptParts,"]", sep = "")))
  temp = temp[!is.na(temp)]
  write("#!/bin/bash", paste("scriptPart", x, ".sh", sep = ""))
  write(temp, paste("scriptPart", x, ".sh", sep = ""), append = T)
  write(paste("bash scriptPart", x, ".sh", " & ", sep = ""), "scriptOfScripts.sh", append = T)
  x = x + 1
}
dir.create("mcast out")
setwd("..")
getwd()
#script of scripts wiggle
setwd(path)
setwd("output")


setwd(path)
setwd("scripts")
getwd()
splitScripts = NULL
scriptStr = unlist(strsplit(wiggleScriptStr, "\n"))
x = 0
scriptParts = ceiling(length(scriptStr)/(threads)) 
write("#!/bin/bash\n", "wiggleScriptOfScripts.sh")
while (x < (threads)){
  temp = eval(parse(text=paste("scriptStr[",x*scriptParts,":",(x+1)*scriptParts,"]", sep = "")))
  temp = temp[!is.na(temp)]
  write("#!/bin/bash", paste("wiggleScriptPart", x, ".sh", sep = ""))
  write(temp, paste("wiggleScriptPart", x, ".sh", sep = ""), append = T)
  write(paste("bash wiggleScriptPart", x, ".sh", " & ", sep = ""), "wiggleScriptOfScripts.sh", append = T)
  x = x + 1
}

#building the wiggle file from wiggleTable
wiggleTable = fread("wiggleTable.txt", header = F)
x = 1
for (i in wiggleTable$V4){
  wiggleTable[x,4] = as.numeric(i)
  x = x + 1
}
x = 1
for (i in wiggleTable$V5){
  wiggleTable[x,5] = as.numeric(i)
  x = x + 1
}
wiggleTable = subset(wiggleTable, !(is.na(wiggleTable$V4)))
wiggleTable = subset(wiggleTable, !(is.na(wiggleTable$V5)))
genes = unique(wiggleTable$V1)
wiggleCoords = NULL
for (i in genes){
  temp = unique(subset(wiggleTable, wiggleTable$V1 == i)$V2)
  start = min(as.numeric(subset(wiggleTable, wiggleTable$V1 == i)$V4))
  end = max(as.numeric(subset(wiggleTable, wiggleTable$V1 == i)$V5))
  sqLength = end - start
  wiggleCoords = rbind(wiggleCoords, c(i, temp[1], start, end, sqLength))
}
wiggleCoords = wiggleCoords[order(wiggleCoords[,2], wiggleCoords[,1]),]

remove(genome) #removing genome from memory to conserve RAM
setwd(path)
setwd("base data")
setwd("phyloP100way")
setwd("chromosomes")
wiggleFile = NULL
for (i in unique(wiggleCoords[,2])){
  chrWiggle = fread(i, header = T)
  temp = subset(wiggleCoords, wiggleCoords[,2] == i)
  for (i in temp[,3]){
    stop = subset(temp, temp[,3] == i)
    #loading chromosome specific wiggle files only temporarily to conserve RAM
    ## wiggle files by UCSC are 1-indexed
    stop = as.numeric(stop[,4]) + 1
    start = as.numeric(i) + 1
    wiggleFile <<- rbind(wiggleFile, chrWiggle[i:stop], use.names = F, col.names = F, row.names = F)
  }
}
setwd(path)
setwd("output")
write.table(wiggleFile, "wiggleTest.wig", col.names= F, row.names = F)

###No longer needed: Change seq files, so that FASTA-header matches the wiggle-coordinates
###probably by first only constructing the coordinate system in the functions above, translating the coordinates to the wiggle coordinates after building those and writing the seq-files as a last step