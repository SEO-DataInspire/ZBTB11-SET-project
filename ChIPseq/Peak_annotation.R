### ChIP-seq peak annotation
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb38<-TxDb.Hsapiens.UCSC.hg38.knownGene

peak.zb<-readPeakFile("ZBTB11_peaks.narrowPeak"); head(peak.zb)
peak.ig<-readPeakFile("IgG_peaks.narrowPeak"); head(peak.ig)
# Annontation
peakAnno.zb<-annotatePeak(peak.zb, tssRegion = c(-3000, 3000),                                 
                          TxDb = txdb38,
                          addFlankGeneInfo=TRUE,flankDistance=5000,
                          annoDb="org.Hs.eg.db")
peakAnno.ig<-annotatePeak(peak.ig, tssRegion = c(-3000, 3000),     
                          TxDb = txdb38,
                          addFlankGeneInfo=TRUE,flankDistance=5000,
                          annoDb="org.Hs.eg.db")
  plotAnnoPie(peakAnno.zb)
  plotAnnoPie(peakAnno.ig)

## output result
write.csv(as.data.frame(peakAnno.zb), file=paste0(gsefilepath,"peakAnno_aZBTB11.csv")) 
write.csv(as.data.frame(peakAnno.ig), file=paste0(gsefilepath,"peakAnno_aIgG.csv")) 


  