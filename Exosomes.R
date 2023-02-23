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
scaleDat<-senscale
senscale <- senscale %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "treatment") %>% 
  pivot_longer(`VIM|P08670`:`TP53BP1|Q12888`,names_to = "protein", values_to = "abundance")


irsearch <- filter(exoDat, grepl(paste(ir$protein, collapse="|"), protein))
rassearch <- filter(exoDat, grepl(paste(ras$protein, collapse="|"), protein))
atvsearch <- filter(exoDat, grepl(paste(atv$protein, collapse="|"), protein))

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

senraw1 <- as.data.frame(senraw) %>% rownames_to_column(var = "prot")


scalemean <- senscale %>% 
  pivot_wider(names_from = treatment, values_from = abundance) %>% 
  column_to_rownames(var = "protein") %>% 
  select(1:16) %>% 
  mutate(mean_ab = rowMeans(., na.rm=TRUE)) %>% 
  rownames_to_column(var = "protein")

proteinGroups<-read.delim(file="chad_16sample_multipdResult_Proteins.txt",sep="\t",header=TRUE) 

########################################################################################################
## ANOVA / DiffEx -- 3 user functions:
##
##   parANOVA.dex      -- create ANOVAout dataframe of tests for differential expression/abundance
##   plotVolc          -- create PDF and HTML Volcano Plots, output volcano settings to variables used later
##   DEXpercentStacked -- create PDF 
##
## By Eric Dammer, Duc Duong, and Qiudong Deng
#########################################
## NOTES
##
## - Output ANOVA+Tukey pairwise stats for volcano and downstream analyses
## - if only 2 comparison groups are present, ANOVA overall p value is equivalent to T test;
##   FDR correction for all proteinwide comparisons provided in that special case (default twoGroupCorrMethod="BH")
## - can be repeat for different subgroup comparisons, if necessary
## - parallelized function. Requires R packages doParallel, parallel, and dependencies
## - writes DEX table typically to pipeline variable ANOVAout and .csv table
## - 0 Tukey values are inaccurate when <1e-9 or -10; Estimates of very small values become 0 in R base Tukey post hoc calculations.
## - So, we implement an option to fallback from Tukey p values less than 1e-10 to Bonferroni-corrected T test (unequal variance)
## - avoid dashes ("-") in strings representing comparison groups (Grouping vector); they will be substituted with '.'
##
##
#########################################
## Required Loaded Data and Parameters ##
#########################################

#ootdir="e:/SysBioPipeline/"
#setwd(rootdir)
#load("SeyfriedPipelineOutput-test.RData")


# Parameters (function may fallback to defaults if unspecified)
Grouping=c("fib", "fib", "fib", "fib", "fib", "fib", "fib", "fib", "fib", "fib", "fib", "fib", "exo", "fib", "exo", "exo")                 # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.
parallelThreads=3                              # number of CPU threads to speed calculation (recommended, 2 or more):
#NETcolors=net$colors                           # list net with slot/vector containing module color assignments; length of vector must be equal to number of rows in cleanDat.
twoGroupCorrMethod="BH"                        # default method for full FDR correction when only 2 groups present is Benjamini-Hochberg; see p.adjust(..., methods= ) options.
outputCSV=TRUE                                 # Output Full Table of statistics?  TRUE/FALSE
outFilePrefix=""                              # typically "4", or step # in pipeline being run; for output file sorting by filename.
outFileSuffix="fib.exo"                 # A description of the project, used as a filename suffix
fallbackIfSmallTukeyP=TRUE                     # Inaccurate Tukey p values < 1e-10 will not replaced with reliable Bonferroni FDR from T test for the pairwise comparison.


source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required input.
#...Tukey p<1e-10 Fallback calculations using Bonferroni corrected T test: 8167 [15%]


## Before plotting volcanoes, choose comparison (column #s) to flip denominator in log2(ratios);
## - any column numbers specified flip the volcano X axis for that pairwise comparison.
head(ANOVAout)
#                F-Value       Pr(>F) copper-control peroxide-control peroxide-copper   diff copper-control diff peroxide-control diff peroxide-copper NETcolors
# VIM|P08670    6.077121 0.0150358128   0.0138057138       0.03492769     0.804518544             1.3474172             1.1415154           -0.2059018 turquoise
# AHNAK|Q09666  2.359204 0.1367470586   0.1230929749       0.24417379     0.855114462             0.8222972             0.6540177           -0.1682795 turquoise
# MYH9|P35579   1.174772 0.3420259215   0.3213307792       0.70687201     0.673431984             0.5022160             0.2680584           -0.2341577     green
# ACTG1|P63261 13.218898 0.0009258414   0.0035954244       0.87534478     0.002040887             0.5481106             0.0653748           -0.4827358 lightcyan
# ACTB|P60709   4.201014 0.0414044802   0.0334517317       0.20512530     0.412588455             0.7005944             0.4401921           -0.2604023 turquoise
# PLEC|Q15149  12.405769 0.0012000067   0.0008621121       0.01779653     0.126763180             0.7427401             0.4842642           -0.2584759    purple
#-------------- ---------- ------------ THESE COLUMNS ARE TUKEY POST-HOC P VALUE COLUMNS FOR PAIRWISE COMPARISONS---  THESE COLUMNS ARE LOG2 MEAN GROUP DIFFERENCES FOR THE SAME PAIRWISE COMPARISONS------

########################################
## Volcano Plots: PDF and HTML        ##
########################################
## Parameters set as variables here.  ##
########################################
## Most will be set to default values ##
## if not explicitly specified.       ##
########################################

FCmin=0                    # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
selectComps=c(3)          # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
flip=c(3)              # ANOVAout column index numbers for p values in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
# As a general rule, the group with less severe effects is usually set to be the denominator (represented by what is right of '-' in ANOVAout column names)
signifP=0.05               # p value threshold for counting Differential Expression points
useNETcolors=FALSE          # use module colors saved to ANOVAout, if available; otherwise when FALSE, specify downColor upColor, and NCcolor (must be valid R color specifications in quotes)
downColor="royalblue"      # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
  upColor="red"              # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
    NCcolor="grey"             # points not significant are this color if useNETcolors=FALSE
      splitColors=FALSE          # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
      highlightGeneProducts=c("CANX|P27824")  # c("APP|P05067","MAPT|P10636","APOE|P02649") ; a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
      labelHighlighted=FALSE     # if true, highlighted spots get text labels with their rownames from ANOVAout
      symbolsOnly=FALSE          # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
      labelTop=0                 # maximum p below which to label all points in the PDF output; OR an integer number of top ranked most significant points to label
      labelSize=4.5              # text label font size, if any labels are found (when labelHighlighted=TRUE or labelTop>0)
      sameScale=FALSE            # When multiple plots are drawn, should they all be the same scale with min and max x and y ranges?
      HTMLout=TRUE               # output interactive HTML copies that can be opened in browser. Requires plotly package.
      outFilePrefix=""          # typically the step # in the pipeline being run
      outFileSuffix="fib.exo"
      # A description of the project, used as a filename suffix
      outputfigs="~/R/Senescence/outputfigs/"        # Location to save figure file output(s)
      
      plotVolc()                 # runs on ANOVAout as input (need not be specified).
      
      ## not run: highlight gene products of interest:
      # highlightGeneProducts=c("APP","SMOC1","MAPT")
      # symbolsOnly=TRUE
      # plotVolc(ANOVAout)
      
      ## not run:  highlight gene products of interest by symbol, and do not use WGCNA module colors:
      # flip=c(3,4,5)
      # useNETcolors=FALSE
      # highlightGeneProducts=c("APP","SMOC1","MAPT")
      # symbolsOnly=TRUE
      # plotVolc(ANOVAout)
      
      
      
      #################################################################
      ## DEx Stacked Bar Plots: PDF(s)                               ##
      #################################################################
      ## Parameters set as variables exported by plotVolc() before.  ##
      #################################################################
      ## Use this function if you have modules in memory created     ##
      ## using the Seyfried Systems Biology Pipeline, after running  ##
      ## parANOVA.dex() and plotVolc() functions.                    ##
      #################################################################
      
      DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.
      
      ###################################################################################################################################
      ## One-Step GSA FET (R piano implementation) WITH USER PARAMETERS
      ##  - by Eric Dammer, Divya Nandakumar
      ## Nicholas Seyfried Lab Bioinformatics - for the lab - 07/27/2022 version 1.01
      ###################################################################################################################################
      ## Preload WGCNA standard pipeline R ression's data to memory, if desired (example minimal RData given here)
      ## otherwise below code runs as a step late in the Seyfried systems biology pipeline,
      ## or following the pipeline's volcano code block, or with .csv input formatted as simple lists
      ## in columns, or the .csv may be kME table output from the global network plots pipeline.
      ###############################################################################################
      
      dev.off()
      
diffex <- read.csv("ANOVA_diffEx-ALL-fib.exo.csv")
diffexo <- diffex %>% 
  filter(diff.fib.exo < 0) %>% 
  filter(fib.exo < 0.05) %>% 
  select(c(X,fib.exo,diff.fib.exo))

diffatlassearch <- filter(diffexo, grepl(paste(atlas$BJ, collapse="|"), X))
diffexomarksearch <- filter(diffexo, grepl(paste(exomarks$Gene, collapse="|"), X))
tempexomarksearch <- filter(rownames_to_column(temp2, var = "X"), grepl(paste(exomarks$Gene, collapse="|"), X))

exoall <- read.delim("EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt", sep = "\t", header = TRUE) %>% 
  filter(SPECIES == "Homo sapiens",
         CONTENT.TYPE == "protein",
         METHODS == "Mass spectrometry")
exoall <- exoall[!duplicated(exoall$GENE.SYMBOL), ]

diffexoallsearch <- filter(exoall, grepl(paste(diffexo$X, collapse="|"), GENE.SYMBOL))

install.packages("ggVennDiagram")
library(ggVennDiagram)

exovenn <- list(Fib = ,
                Exo = )
