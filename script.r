library('data.table')
#used for "fread"
library("liftOver")
#library('ggplot2')
path = getwd()
path = paste(path,"/Documents/Pankreaspraktikum/data", sep="")
setwd(path)
#testPath = paste(getwd(),"ABC/AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv", sep="")
dataABC = fread("ABC/AllPredictions.AvgHiC.ABC0.02.ModelRegions.csv")
#Datenbankoutput der ABC Max Methode
genes = fread("genesWitt.txt", header = FALSE)
#Vorl?ufige Liste an Keywords (Gene und Loci), zusammengeschrieben aus der Publikationsliste von Prof Witt 
preFilteredGenes = subset(dataABC,dataABC$TargetGene %in% genes$V1)
filteredGenes = subset(preFilteredGenes,preFilteredGenes$CellType == "body_of_pancreas-ENCODE")
foundGenes = subset(genes, !(genes$V1 %in% filteredGenes$TargetGene))
#Filterung der Datenbank auf Eintraege in denen "TargetGenes" auf der Liste der relevanten stehen
#View(filteredGenes)
TFBSs = fread("Human_TF_MotifList_v_1.01.csv")
getexMeans = fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")
#getex = fread(paste(path, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep=""))
getex = getexMeans
tfTPM = subset(getex, getex$Description %in% TFBSs$`HGNC symbol`)
#plot(tfTPM$Pancreas)
tfTranscribedInP = subset(tfTPM, tfTPM$Pancreas > 1)
tfFilterListe = fread("TFsWitt.txt", header = F)
#hgncToCisbp = subset(TFBSs, TFBSs$`HGNC symbol` %in% diseasesPancreatitisSignificantTFs$V2)
hgncToCisbp = subset(TFBSs, TFBSs$`HGNC symbol` %in% tfFilterListe$V1)
hgncToCisbpList = subset(hgncToCisbp, hgncToCisbp$`HGNC symbol` %in% tfTranscribedInP$Description)$"CIS-BP ID"

diseases = fread("human_disease_textmining_full.tsv")
diseasesPancreatitis = subset(diseases, diseases$V3 == "DOID:4989")
diseasesPancreatitisSignificant = subset(diseasesPancreatitis, diseasesPancreatitis$V5 > 2)
diseasesPancreatitisSignificantTFs = subset(diseasesPancreatitisSignificant, diseasesPancreatitisSignificant$V2 %in% tfTranscribedInP$Description)
write(diseasesPancreatitisSignificantTFs$V2, file = "pancreatitisTFs.txt")


#conversion script fuer individuelle Motifs
# setwd("PWMs/")
# setwd("txt files/")
# fileList = list.files()
# setwd(gsub("/txt files","",getwd()))
# 
# for (i in fileList){
#   setwd("txt files/")
#   fileName = i
#   fileName = gsub(".txt","", fileName)
#   motif = fread(i)
#   pwm=motif[,2:5]
#   pwmNameless=unname(pwm)
#   identifier = fileName
#   alternateName =paste(fileName, "testest", sep="")
#   setwd(gsub("/txt files","",getwd()))
#   write("MEME version 4", file=paste(fileName, ".meme", sep=""))
#   write("ALPHABET= ACGT", file=paste(fileName, ".meme", sep=""), append =T)
#   write(paste("MOTIF", identifier, alternateName, sep=" "), file=paste(fileName, ".meme", sep=""), append=T)
#   write("letter-probability matrix:", file=paste(fileName, ".meme", sep=""), append=T)
#   write.table(pwm,file=paste(fileName, ".meme", sep=""), append= T, col.names= F, row.names = F)
# }

#conversion script into one file
setwd("PWMs/")
setwd("txt files/")
fileList = list.files()
fileList = subset(fileList, gsub(".txt","",fileList) %in% hgncToCisbpList)
setwd(gsub("/txt files","",getwd()))
write("MEME version 4", file=paste("master.meme", sep=""))
write("ALPHABET= ACGT", file=paste("master.meme", sep=""), append =T)

for (i in fileList){
  #TODO: alternateName sollte HGNC symbol für korrespondierendes Gen sein
  #TODO: pruefen ob das jeweilige motif die bedingungen von SEA an eine PSPM erfuellt (keine probability = 1, keine summe != 1)
  setwd("txt files/")
  fileName = "master"
  motifName = i
  motifName = gsub(".txt","", motifName)
  motif = fread(i)
  pwm=motif[,2:5]
  pwmNameless=unname(pwm)
  identifier = motifName
  alternateName = subset(hgncToCisbp, hgncToCisbp$`CIS-BP ID` == gsub(".txt","",i))$"HGNC symbol"[1]
  print(alternateName)
  setwd(gsub("/txt files","",getwd()))

  write(paste("MOTIF", identifier, alternateName, sep=" "), file=paste(fileName, ".meme", sep=""), append=T)
  write("letter-probability matrix:", file=paste(fileName, ".meme", sep=""), append=T)
  write.table(pwm,file=paste(fileName, ".meme", sep=""), append= T, col.names= F, row.names = F)
  write("\n", file=paste(fileName, ".meme", sep=""), append=T)
}
setwd(gsub("/PWMs","",getwd()))

setwd("hg38/FASTA")
chromosomes = paste("chr", c(1:22), sep="") 
chromosomes = c(chromosomes, "chrX", "chrY")
genome = NULL
for (i in chromosomes){
  genome[i] = fread(paste(i,".fna", sep=""))
}  

setwd(gsub("/hg38/FASTA","", getwd()))

#liftover

shiftedGenes = filteredGenes
for (i in filteredGenes$start){
  filteredGenes$start = filteredGenes$start -1 #enumeration shift
  filteredGenes$end = filteredGenes$end -1
}

newFormat = makeGRangesFromDataFrame(shiftedGenes)
liftedSeqs = liftOver(newFormat, import.chain(system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")))

setwd(gsub("/hg38/FASTA","",getwd()))
setwd("sequences")

#script für ein sammelfile fuer alle seqs
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
  write(paste(">",start," Homo sapiens ", subset(filteredGenes, filteredGenes$name ==i)$chr,", GRCh38.p14", sep=""), file="master.fa", append = T)
  write(targetSequence, file="master.fa", append = T)
  write("", file="master.fa", append = T)
  x = x+1
  print(x)

}
setwd(gsub("/sequences","",getwd()))

#script fuer individuelle files je sequenz
#x = 1
# for (i in filteredGenes$name){
#   start = start(liftedSeqs[[x]])
#   end = end(liftedSeqs[[x]])
#   startRow = floor(start/80)+1
#   startPoint = start - startRow*80 +80
#   endRow = floor(end/80)+1
#   endPoint = end - endRow*80 +80
#   targetRows = paste(unlist(eval(parse(text = paste("genome$",subset(filteredGenes, filteredGenes$name ==i)$chr,"[",startRow,":",endRow,"]", sep="")))), collapse="")
#   targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
#   write(paste(">",start," Homo sapiens ", subset(filteredGenes, filteredGenes$name ==i)$chr,", GRCh38.p14", sep=""), file=paste(i, ".fa", sep=""))
#   write(targetSequence, file=paste(i, ".fa", sep=""), append = T)
#   x = x+1
#   print(x)
#   
# }
# setwd(gsub("/sequences","",getwd()))
# setwd(gsub("/hg38/FASTA","",getwd()))
# setwd("sequences")
# 
# for (i in filteredGenes$name){
#   start = subset(filteredGenes, filteredGenes$name == i)$start
#   end = subset(filteredGenes, filteredGenes$name == i)$end
#   startRow = floor(start/80)+1
#   startPoint = start - startRow*80 +80
#   endRow = floor(end/80)+1
#   endPoint = end - endRow*80 +80
#   targetRows = paste(unlist(eval(parse(text = paste("genome$",subset(filteredGenes, filteredGenes$name ==i)$chr,"[",startRow,":",endRow,"]", sep="")))), collapse="")
#   targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))
#   write(paste(">",start," Homo sapiens ", subset(filteredGenes, filteredGenes$name ==i)$chr,", GRCh38.p14", sep=""), file=paste(i, ".fa", sep=""))
#   write(targetSequence, file=paste(i, ".fa", sep=""), append = T)
#   
# }
# start = 115913256
# end = 115914368
# startRow = floor(start/80)+1
# startPoint = start - startRow*80 +80
# endRow = floor(end/80)+1
# endPoint = end - endRow*80 +80
# targetRows = paste(unlist(genome$chr7[startRow:endRow]), collapse="")
# targetSequence = substr(targetRows, startPoint, (nchar(targetRows)-(80-endPoint)))



