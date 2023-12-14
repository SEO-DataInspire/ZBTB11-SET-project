### RNA-seq analysis

library(DESeq2)
library(dplyr)

siCtr.rep1<-read.table('si-Ctr_pcgene_rep1.txt', header = T, row.names = 1,stringsAsFactors = F)  
siCtr.rep2<-read.table('si-Ctr_pcgene_rep2.txt', header = T, row.names = 1,stringsAsFactors = F) 
siZBTB11.rep1<-read.table('si-ZBTB11_pcgene_rep1.txt', header = T, row.names = 1,stringsAsFactors = F) 
siZBTB11.rep1<-read.table('si-ZBTB11_pcgene_rep2.txt', header = T, row.names = 1,stringsAsFactors = F) 
siSET.rep1<-read.table('si-SET_pcgene_rep1.txt', header = T, row.names = 1,stringsAsFactors = F) 
siSET.rep1<-read.table('si-SET_pcgene_rep2.txt', header = T, row.names = 1,stringsAsFactors = F) 

#--------------- DEG analysis of ZBTB11 --------------#

## Count data preparation
 data.ZBTB11<-data.frame(siCtr_1=siCtr.rep1[,1], siCtr_2=siCtr.rep2[,1],
 						 siZBTB11_1=siZBTB11.rep1[,1],siZBTB11_2=siZBTB11.rep2[,1],
						 row.names=rownames(siCtr_1))

 data.ZBTB11<-data.ZBTB11[rowSums(data.ZBTB11==0)==0,] # excloud genes which value is zero 
  
 SampleGroup<-data.frame(condition=c("NC","NC","KD","KD"),                          
  	   					 treat=c("NC","NC", "KD","KD"))
 cts<-data.ZBTB11
 colData<-SampleGroup
 dds<-DESeqDataSetFromMatrix(countData = cts,
                             colData = colData,
                             design = ~condition)
dds
dds$condition<-factor(dds$condition, levels = c("NC", "KD"))

## Differential expression analysis by count data
  set.seed(123)
  dds<-DESeq(dds)
  res<-results(dds)
  head(res)

## output result
  write.csv(res, "DEG_analysis_of_siZBTB11.csv")

#--------------- DEG analysis of SET --------------#

## Count data preparation
 data.SET<-data.frame(siCtr_1=siCtr.rep1[,1], siCtr_2=siCtr.rep2[,1],
 						 siSET_1=siSET.rep1[,1],siSET_2=siSET.rep2[,1],
						 row.names=rownames(siCtr_1))

 data.SET<-data.SET[rowSums(data.SET==0)==0,] # excloud genes which value is zero 
  
 SampleGroup<-data.frame(condition=c("NC","NC","KD","KD"),                          
  	   					 treat=c("NC","NC", "KD","KD"))
 cts<-data.SET
 colData<-SampleGroup
 dds<-DESeqDataSetFromMatrix(countData = cts,
                             colData = colData,
                             design = ~condition)
dds
dds$condition<-factor(dds$condition, levels = c("NC", "KD"))

## Differential expression analysis by count data
  set.seed(123)
  dds<-DESeq(dds)
  res<-results(dds)
  head(res)

## output result
  write.csv(res, "DEG_analysis_of_siSET.csv")
 



