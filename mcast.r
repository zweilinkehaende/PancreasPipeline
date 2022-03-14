library("data.table")
library("liftOver")
# library("biomaRt")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("org.Hs.eg.db")

test <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
getGene("HNF1A", type = 'hgnc_symbol', mart)
setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()
setwd("base data")
# HGNCtoENS = fread("ens_hgnc")

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
hg38Annotated = (TxDb.Hsapiens.UCSC.hg38.knownGene)
#defines function to create seq file
geneGroup = "A.gene"
kbRange = 10000
prepMcast = function(geneGroup, kbRange){
  #requires libraries "data.table" and "liftOver"
  #requires global variables "dataABC", "tfTranscribedInP" and "genome"
  
  setwd("input")
  geneListOfGroup = fread(geneGroup, header = FALSE)
  #gene deren sequenzen zu untersuchen sind

  setwd(gsub("/input","", getwd()))

  fileName = gsub(".gene", "", genes)
  fileName = paste(fileName, ".fa", sep="")
  ensList = NULL
  for (i in geneListOfGroup$V1){
    temp = mapIds(org.Hs.eg.db, keys = i,column = "ENTREZID", keytype = "SYMBOL")
    ensList = rbind(ensList, temp[[1]])
  }
  test = select(TxDb.Hsapiens.UCSC.hg38.knownGene, column = "GENEID", keys = ensList, keytype = "GENEID")
  
  #script fuer ein sammelfile fuer alle seqs
  x = 1
  for (i in ensList){
    temp = subset(hg38Annotated, hg38Annotated$ == i)
    tss = eval(parse(text = paste("start(",subset(hg38Annotated, hg38Annotated$tx_name == i),"[i))", sep = "")))
    start = tss - kbLength
    end = tss + kbLength
    startRow = floor(start/80)+1 #+1 wegen floor()
    startPoint = start - startRow*80 + 80
    endRow = floor(end/80)+1
    endPoint = end - endRow*80 + 80
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
    identifier = motifName
    alternateName = subset(hgncToCisbp, hgncToCisbp$`CIS-BP ID` == gsub(".txt","",j))$"HGNC symbol"[1]
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
regionsList = c(1:10)
for (i in geneList){
  for (j in regionsList){
    prepMcast(i, (j*10000))
  }
  
}

##writes the script to run sea for all combinations
setwd("output")

library(parallel)
threads = detectCores() -1 #using 1 thread fewer than available to prevent the computer from freezing up
dirs = subset(list.files(), !(grepl(".meme", list.files())))
scriptStr = NULL
for (i in dirs){
  sharedFiles = subset(list.files(), grepl(".meme", list.files()))
  
  setwd(paste(i))
  files = list.files()
  scriptStr = paste(scriptStr, 'sea --p "', i, "/", subset(files, (grepl(".fa", files))), '" --m "', subset(files, (grepl(".meme", files))), '"\n', sep = "")
  setwd("..")
}
setwd("..")
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
