library(tidyverse)
library(readxl)
load("SenescenceProteomicsCleanNumeric.RData")
ir <- read_xlsx("Supplemental_Table_1.xlsx", sheet = "IR SASP", cell_cols("B"))
ras<- read_xlsx("Supplemental_Table_1.xlsx", sheet = "RAS SASP", cell_cols("B"))
atv<- read_xlsx("Supplemental_Table_1.xlsx", sheet = "ATV SASP", cell_cols("B"))
exomarks <- read.table("exomarks.txt", fill = TRUE, header = TRUE) %>% 
  select(Gene)
bca <- read_xlsx("Hales16lysates_120522.xlsx", sheet = "Processing", cell_rows(c(2:18))) %>% 
  select("Sample number":"Protein conc (µg/µL)") %>% 
  mutate(scale = 2.654159/`Amount digested  (µg)`)
senraw <- read.csv("sen.PD.LFQ_rawAbundances.NoCorrection.csv") %>% 
  column_to_rownames(var = "X") %>% 
  as.matrix()
colnames(senraw) = 1:16
ir<-strsplit(ir$UniProtIds, ";",) %>% 
  unlist() %>% 
  as.data.frame()
names(ir) <- "protein"
ras<-strsplit(ras$UniProtIds, ";",) %>% 
  unlist() %>% 
  as.data.frame()
names(ras) <- "protein"
atv<-strsplit(atv$UniProtIds, ";",) %>% 
  unlist() %>% 
  as.data.frame()
names(atv) <- "protein"
exoDat <- as.data.frame(subset(cleanDat, select = c(`Abundance..F13`,`Abundance..F15`,`Abundance..F16`))) %>% 
  rename(Copper=`Abundance..F13`,
         Peroxide=`Abundance..F15`,
         Control=`Abundance..F16`)
LThalfSamples<-length(colnames(exoDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(exoDat)) %% 2)==1) { 0.5 } else { 1.0 }
temp2<-as.data.frame(exoDat[which(rowSums(as.matrix(is.na(exoDat)))>LThalfSamples),])
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(exoDat)[which(rowSums(as.matrix(is.na(exoDat)))>LThalfSamples)]
}
if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(exoDat)))>LThalfSamples);rowsRemoved<-rownames(exoDat)[IndexHighMissing]; exoDat<-exoDat[-IndexHighMissing,]; }
exoDat <- exoDat %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "treatment") %>% 
  pivot_longer(`VIM|P08670`:`TP53BP1|Q12888`,names_to = "protein", values_to = "abundance")
fibDat <- as.data.frame(subset(cleanDat, select = c(`Abundance..F1`:`Abundance..F12`,`Abundance..F14`))) %>% 
  mutate(`Abundance..F13`=`Abundance..F14`, .before = `Abundance..F14`) %>% 
  mutate(`Abundance..F15`=`Abundance..F14`)
LThalfSamples<-length(colnames(fibDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(fibDat)) %% 2)==1) { 0.5 } else { 1.0 }
temp2<-as.data.frame(fibDat[which(rowSums(as.matrix(is.na(fibDat)))>LThalfSamples),])
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(fibDat)[which(rowSums(as.matrix(is.na(fibDat)))>LThalfSamples)]
}
if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(fibDat)))>LThalfSamples);rowsRemoved<-rownames(fibDat)[IndexHighMissing]; fibDat<-fibDat[-IndexHighMissing,]; }
fibDat <- fibDat %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "treatment") %>% 
  pivot_longer(`VIM|P08670`:`LRRC40|Q9H9A6`,names_to = "protein", values_to = "abundance")
senraw[is.na(senraw)] <- 0
senscale <- senraw%*%diag(bca$scale)
senscale[senscale == 0] <- NA
senscale <- senscale %>% 
  log(base = 2)
colnames(senscale) = 1:16
LThalfSamples<-length(colnames(senscale))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(senscale)) %% 2)==1) { 0.5 } else { 1.0 }
temp2<-as.data.frame(senscale[which(rowSums(as.matrix(is.na(senscale)))>LThalfSamples),])
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(senscale)[which(rowSums(as.matrix(is.na(senscale)))>LThalfSamples)]
}
if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(senscale)))>LThalfSamples);rowsRemoved<-rownames(senscale)[IndexHighMissing]; senscale<-senscale[-IndexHighMissing,]; }
senscale <- senscale %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "treatment") %>% 
  pivot_longer(`VIM|P08670`:`TP53BP1|Q12888`,names_to = "protein", values_to = "abundance")


irsearch <- filter(exoDat, grepl(paste(ir$protein, collapse="|"), protein))
rassearch <- filter(exoDat, grepl(paste(ras$protein, collapse="|"), protein))
atvsearch <- filter(exoDat, grepl(paste(atv$protein, collapse="|"), protein))

exoDat %>% 
  filter(str_detect(protein, "Q9H772"))
fibDat %>% 
  filter(str_detect(protein, "Q9H772"))
#GREM2


#Exosome Enrichment
#16 markers searched
#8 found
#1284 total proteins
exoDat %>% 
  filter(str_detect(protein, "P08962"))
fibDat %>% 
  filter(str_detect(protein, "P08962"))
#CD63
exoDat %>% 
  filter(str_detect(protein, "P60033"))
fibDat %>% 
  filter(str_detect(protein, "P60033"))
#CD81
exoDat %>% 
  filter(str_detect(protein, "Q8WUM4"))
fibDat %>% 
  filter(str_detect(protein, "Q8WUM4"))
#ALIX
exoDat %>% 
  filter(str_detect(protein, "P04083"))
fibDat %>% 
  filter(str_detect(protein, "P04083"))
#ANXA1
exoDat %>% 
  filter(str_detect(protein, "P07355"))
fibDat %>% 
  filter(str_detect(protein, "P07355"))
#ANXA2
exoDat %>% 
  filter(str_detect(protein, "P08758"))
fibDat %>% 
  filter(str_detect(protein, "P08758"))
#ANXA5
exoDat %>% 
  filter(str_detect(protein, "P02751"))
fibDat %>% 
  filter(str_detect(protein, "P02751"))
#FN1
exoDat %>% 
  filter(str_detect(protein, "HSP90"))
fibDat %>% 
  filter(str_detect(protein, "HSP90"))
#HSP90

exosearch <- filter(exoDat, grepl(paste(exomarks$Gene, collapse="|"), protein)) %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  select(Copper:Control) %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")
fibsearch <- filter(fibDat, grepl(paste(exomarks$Gene, collapse="|"), protein)) %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  select(`Abundance..F1`:`Abundance..F15`) %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")
compsearch <- filter(fibsearch, grepl(paste(exosearch$protein, collapse="|"), protein))
scalesearch <- filter(senscale, grepl(paste(exomarks$Gene, collapse="|"), protein)) %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  select(1:16) %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")


atlas <- read_xlsx("core_atlas.xlsx", sheet = "Supplementary Table 1", cell_cols("C"))
exoatlassearch <- filter(exoDat, grepl(paste(atlas$BJ, collapse="|"), protein)) %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")
fibatlassearch <- filter(fibDat, grepl(paste(atlas$BJ, collapse="|"), protein)) %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")
scaleatlassearch <- filter(senscale, grepl(paste(atlas$BJ, collapse="|"), protein)) %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")
rm(temp2)
rm(IndexHighMissing)
rm(LThalfSamples)
rm(rowsRemoved)
