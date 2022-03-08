#script
library('data.table')
#used for "fread"
library("liftOver")
path = getwd()
setwd("base data")

dataABC = fread("AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv")
preFilteredGenes = subset(dataABC, dataABC$CellType == "body_of_pancreas-ENCODE")
#Datenbankoutput der ABC Max Methode
TFBSs = fread("Human_TF_MotifList_v_1.01.csv")
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

prep = function(relevantGenes, relevantTFs, sqName, PSPMname, sqFilter){
  #requires libraries "data.table" and "liftOver"
  #requires global variables "dataABC", "tfTranscribedInP" and "genome"
  
  # setwd("input")
  setwd("input")
  genes = fread(relevantGenes, header = FALSE)
  #gene deren sequenzen zu untersuchen sind
  filteredGenes = subset(preFilteredGenes,preFilteredGenes$TargetGene %in% genes$V1)
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
    whitelist = rbind(whitelist, filteredGenes = subset(filteredGenes, filteredGenes$class == "promotor"))
  }

  # filteredGenes = subset(filteredGenes, !(filteredGenes$class == notThisSq))


  tfFilterListe = fread(relevantTFs, header = F)
  hgncToCisbp = subset(TFBSs, TFBSs$`HGNC symbol` %in% tfFilterListe$V1)
  hgncToCisbpList = subset(hgncToCisbp, hgncToCisbp$`HGNC symbol` %in% tfTranscribedInP$Description)$"CIS-BP ID"
  
  # setwd(gsub("/input","", getwd()))
  setwd(gsub("/input","", getwd()))

  #conversion script into one file
  setwd("base data/PSPMs")
  fileList = list.files()
  fileList = subset(fileList, gsub(".txt","",fileList) %in% hgncToCisbpList)
  setwd(gsub("base data/PSPMs","", getwd()))
  
  setwd("output")
  outputDirName = paste(gsub(".gene", "", relevantGenes), gsub(".tf", "", relevantTFs), sqFilter)
  dir.create(outputDirName)
  setwd(outputDirName)
  fileName = PSPMname
  write("MEME version 4", file= paste(fileName,".meme", sep=""))
  write("ALPHABET= ACGT", file= paste(fileName,".meme", sep=""), append =T)
  
  for (i in fileList){
    #TODO: pruefen ob das jeweilige motif die bedingungen von SEA an eine PSPM erfuellt (keine probability = 1, keine summe != 1)
    setwd(gsub(paste("/", outputDirName, sep=""),"", getwd()))
    setwd(gsub("/output","",getwd()))
    setwd("base data/PSPMs")
    motifName = i
    motifName = gsub(".txt","", motifName)
    motif = fread(i)
    pwm = motif[,2:5]
    pwmNameless = unname(pwm)
    identifier = motifName
    alternateName = subset(hgncToCisbp, hgncToCisbp$`CIS-BP ID` == gsub(".txt","",i))$"HGNC symbol"[1]
    setwd(gsub("base data/PSPMs","",getwd()))
    setwd("output")
    setwd(outputDirName)
    
    write(paste("MOTIF", identifier, alternateName, sep=" "), file=paste(fileName, ".meme", sep=""), append=T)
    write("letter-probability matrix:", file=paste(fileName, ".meme", sep=""), append=T)
    write.table(pwm,file=paste(fileName, ".meme", sep=""), append= T, col.names= F, row.names = F)
    write("\n", file=paste(fileName, ".meme", sep=""), append=T)
  }

  #liftover
  
  shiftedGenes = filteredGenes
  for (i in filteredGenes$start){
    filteredGenes$start = filteredGenes$start -1 #enumeration shift
    filteredGenes$end = filteredGenes$end -1
  }
  
  newFormat = makeGRangesFromDataFrame(shiftedGenes)
  liftedSeqs = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))
  
  fileName = sqName
  fileName = paste(fileName, ".fa", sep="")
  #script fuer ein sammelfile fuer alle seqs
  x = 1
  for (i in filteredGenes$name){
    start = start(liftedSeqs[[x]])
    end = end(liftedSeqs[[x]])
    startRow = floor(start/80)+1 #+1 wegen floor()
    startPoint = start - startRow*80 +80
    endRow = floor(end/80)+1
    endPoint = end - endRow*80 +80
    targetRows = paste(unlist(eval(parse(text = paste("genome$",subset(filteredGenes, filteredGenes$name ==i)$chr,"[",startRow,":",endRow,"]", sep="")))), collapse="")
    targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
    write(paste(">",start," Homo sapiens ", subset(filteredGenes, filteredGenes$name ==i)$chr,", GRCh38.p14", sep=""), file= fileName, append = T)
    write(targetSequence, file= fileName, append = T)
    write("", file=fileName, append = T)
    x = x+1
    
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


for (i in geneList){
  for (j in tfList){
    prep(i, j, gsub(".gene", "", i), gsub(".tf", "", j), "a")
  }
}
setwd("output")
dirs = list.files()
scriptStr = NULL
for (i in dirs){
  scriptStr = paste(scriptStr, paste('cd "', i, '"\n', sep = ""), sep = "")
  setwd(paste(i))
  files = list.files()
  scriptStr = paste(scriptStr, paste('sea --p "', subset(files, (grepl(".fa", files))), '" --m "', subset(files, (grepl(".meme", files))), '"\n', sep = ""), sep = "")
  scriptStr = paste(scriptStr, "cd ..\n", sep = "")
  setwd("..")
  }
setwd("..")
write(scriptStr, "bashScript")
