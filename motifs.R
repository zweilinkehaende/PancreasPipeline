############
#builds the shared .meme motif files
library("data.table")
setwd("Documents/Pankreaspraktikum/script_input_output")
path = getwd()
setwd("base data")

TFBSs = fread("Human_TF_MotifList_v_1.01.csv")
TFBSs = subset(TFBSs, TFBSs$`Best Motif(s)? (Figure 2A)`== "TRUE")
getex = fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")
tfTPM = subset(getex, getex$Description %in% TFBSs$`HGNC symbol`)
tfTranscribedInP = subset(tfTPM, tfTPM$Pancreas > 1)

setwd(path)
setwd("input")
tfList = list.files()
setwd("..")
tfFilterListe = subset(tfList, grepl(".tf", tfList))
for (i in tfFilterListe){
  # i = "allTFs.tf"
  setwd("input")
  temp = fread(i, header = F)
  setwd("..")
  if (!(i == "allTFs.tf")){
    hgncToCisbp = subset(TFBSs, TFBSs$`HGNC symbol` %in% temp$V1)
  }
  if (i == "allTFs.tf"){
    hgncToCisbp = TFBSs
  }
  hgncToCisbpList = subset(hgncToCisbp, hgncToCisbp$`HGNC symbol` %in% tfTranscribedInP$Description)$"CIS-BP ID"
  setwd("base data/PSPMs")
  fileList = list.files()
  fileList = subset(fileList, gsub(".txt","",fileList) %in% hgncToCisbpList)
  setwd("..")
  setwd("..")
  
  setwd("common files")
  fileName = gsub(".tf", "", i)
  write("MEME version 4", file= paste(fileName,".meme", sep=""))
  write("ALPHABET= ACGT", file= paste(fileName,".meme", sep=""), append =T)
  setwd("..")
  
  for (j in fileList){
    #TODO: pruefen ob das jeweilige motif die bedingungen von SEA an eine PSPM erfuellt (keine probability = 1, keine summe != 1)
    setwd("base data/PSPMs")
    # j = "HKR1.RCADE.txt"
    motifName = j
    motifName = gsub(".txt","", motifName)
    motif = fread(j)
    pwm = motif[,2:5]
    pwmNameless = unname(pwm)
    identifier = subset(hgncToCisbp, hgncToCisbp$`CIS-BP ID` == gsub(".txt","",j))$"HGNC symbol"[1]
    alternateName = motifName
    setwd(gsub("base data/PSPMs","",getwd()))
    
    setwd("common files")
    write(paste("MOTIF", identifier, alternateName, sep=" "), file=paste(fileName, ".meme", sep=""), append=T)
    write("letter-probability matrix:", file=paste(fileName, ".meme", sep=""), append=T)
    write.table(pwm,file=paste(fileName, ".meme", sep=""), append= T, col.names= F, row.names = F)
    write("\n", file=paste(fileName, ".meme", sep=""), append=T)
    setwd("..")
  }
}

#TODO: 