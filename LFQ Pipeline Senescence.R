library(tidyverse)
options(stringsAsFactors=FALSE)

rootdir="~/R/Senescence/"
setwd(rootdir)
load("SenescenceProteomicsCleanNumeric.RData")  ##cleanDat has proteins on row and samples on column headers
outputfigs="~/R/Senescence/outputfigs"
outputtabs="~/R/Senescence/outputtabs"

#FibPellets #cheated controls... :[[[[[
cleanDat <- as.data.frame(subset(cleanDat, select = c(`Abundance..F1`:`Abundance..F12`,`Abundance..F14`))) %>% 
  mutate(`Abundance..F13`=`Abundance..F14`, .before = `Abundance..F14`) %>% 
  mutate(`Abundance..F15`=`Abundance..F14`)
numericMeta <- numericMeta %>% filter(type == "fib") %>% 
  mutate(count = c(1,1,1,1,1,1,1,1,1,1,1,1,3)) %>%
  uncount(count) %>% 
  mutate(sample=c(1:15))

library("doParallel")
stopCluster(clusterLocal) #in case already started.
parallelThreads=3

clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)
library(WGCNA)
enableWGCNAThreads() #speeds outlier check calculations

#=============================#
#  Check and Remove Outliers  # 
#=============================#
sdout=3 #Z.k SD fold for outlier threshold
outliers.noOLremoval<-outliers.All<-vector()
cleanDat.noOLremoval<-cleanDat

for (repeated in 1:20) {
  normadj <- (0.5+0.5*bicor(cleanDat,use="pairwise.complete.obs")^2)
  
  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- (ku-mean(ku))/sqrt(var(ku))  #corrected, moved parenthesis open to leading position
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku < mean(z.ku)-sdout*sd(z.ku))  #previously had | z.ku > mean(z.ku)+sdout*sd(z.ku)) to remove outliers on high connectivity end...
  print(paste0("There are ",sum(outliers)," outlier samples based on a bicor distance sample network connectivity standard deviation above ",sdout,".  [Round ",repeated,"]"))
  targets.All=numericMeta
  
  cleanDat <- cleanDat[,!outliers] 
  numericMeta <- targets <- targets.All[!outliers,]
  outliers.All<-c(outliers.All,outliers)
  
  if (sum(outliers)==0) break;
} #repeat up to 20 times

#All outliers removed
print(paste0("There are ",sum(outliers.All)," total outlier samples removed in ",repeated," iterations:"))
names(which(outliers.All))
outliersRemoved<-names(which(outliers.All))
#Note outliers as comment below, copied from R session.
#[1] "There are 0 total outlier samples removed in 1 iterations:"

## Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }

## If operating on log2(FPKM) data, remove rows with >=50% originally 0 FPKM values (only if there are some rows to be removed)
#IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
#temp2<-data.frame(ThrowOut=apply(cleanDat,1,function(x) length(x[x==log2(0+0.05)])>LThalfSamples))
#cleanDat<-cleanDat[!temp2$ThrowOut,]
#dim(cleanDat) #still have x genes, now for y total samples

## If working on log2(protein abundance or ratio) with NA missing values; Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a data frame)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples);rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }
#remove rows with zero variance
adj <- adjacency(t(as.matrix(cleanDat)), type = "signed", power = 12.5, corFnc = "bicor")
adjnona <- as.data.frame(adj) %>% 
  filter(if_any(.cols = everything(), ~ is.na(.))) %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(if_any(.cols = everything(), ~ is.na(.)))
temp3 <- cleanDat %>% 
  filter(row.names(cleanDat) %in% names(adjnona))
if (nrow(temp3)>0) { IndexZeroVar=which(row.names(cleanDat) %in% names(adjnona));rowsRemovedvars<-rownames(cleanDat)[IndexZeroVar]; cleanDat<-cleanDat[-IndexZeroVar,]; }

rowsRemoved
rowsRemovedvars

dim(cleanDat)
#[1] 2001   15


## Estimate SFT power
unloadNamespace("WGCNA")
library(WGCNA)

stopCluster(clusterLocal) #in case already started.
parallelThreads=3

clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)
allowWGCNAThreads() #speeds pickSoftThreshold function and blockwiseModules fn below, for some reason enableWGCNAthreads() does not.

powers <- seq(2,20,by=0.5)
tableSFT<-pickSoftThreshold(t(cleanDat),powerVector=powers,corFnc="bicor",blockSize=50000,verbose=3,networkType="signed") #[[2]]
plot(tableSFT[[2]][,1],tableSFT[[2]][,2],xlab="Power (Beta)",ylab="SFT R?")

# Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
# 1    2.0   0.9460  1.5900         0.9360   979.0    1070.0   1340
# 2    2.5   0.8920  1.1200         0.8610   857.0     943.0   1250
# 3    3.0   0.7830  0.7930         0.7220   759.0     838.0   1170
# 4    3.5   0.5830  0.5270         0.4690   678.0     748.0   1100
# 5    4.0   0.3670  0.3410         0.1920   610.0     672.0   1040
# 6    4.5   0.1610  0.1980        -0.0601   553.0     606.0    986
# 7    5.0   0.0212  0.0698        -0.2340   504.0     547.0    938
# 8    5.5   0.0103 -0.0551        -0.2190   462.0     497.0    896
# 9    6.0   0.0658 -0.1560        -0.1320   425.0     453.0    857
# 10   6.5   0.1280 -0.2480        -0.0304   393.0     412.0    822
# 11   7.0   0.1710 -0.3280         0.0433   365.0     378.0    790
# 12   7.5   0.2160 -0.4070         0.1050   339.0     348.0    759
# 13   8.0   0.2560 -0.4780         0.1600   316.0     319.0    731
# 14   8.5   0.2840 -0.5300         0.2090   296.0     294.0    705
# 15   9.0   0.3150 -0.5640         0.2550   278.0     270.0    681
# 16   9.5   0.3360 -0.5900         0.2850   261.0     249.0    658
# 17  10.0   0.3460 -0.6410         0.3170   246.0     231.0    636
# 18  10.5   0.3640 -0.6710         0.3480   232.0     215.0    616
# 19  11.0   0.3910 -0.6950         0.3840   219.0     200.0    598
# 20  11.5   0.3930 -0.7380         0.4000   208.0     186.0    580
# 21  12.0   0.4170 -0.7650         0.4280   197.0     172.0    563
# 22  12.5   0.4380 -0.7800         0.4590   187.0     161.0    547
# 23  13.0   0.4560 -0.7940         0.4810   178.0     151.0    532
# 24  13.5   0.4730 -0.8110         0.5010   170.0     141.0    518
# 25  14.0   0.4880 -0.8220         0.5190   162.0     132.0    504
# 26  14.5   0.5030 -0.8350         0.5330   154.0     123.0    490
# 27  15.0   0.5020 -0.8630         0.5370   148.0     115.0    478
# 28  15.5   0.5020 -0.8960         0.5480   141.0     108.0    466
# 29  16.0   0.5130 -0.9040         0.5630   135.0     102.0    454
# 30  16.5   0.5250 -0.9140         0.5770   130.0      95.6    443
# 31  17.0   0.5340 -0.9170         0.5910   124.0      89.9    432
# 32  17.5   0.5620 -0.9130         0.6180   119.0      84.5    422
# 33  18.0   0.5760 -0.9290         0.6360   115.0      79.5    412
# 34  18.5   0.5870 -0.9410         0.6480   110.0      74.9    402
# 35  19.0   0.5970 -0.9500         0.6580   106.0      70.6    393
# 36  19.5   0.6060 -0.9520         0.6740   102.0      66.6    384
# 37  20.0   0.6130 -0.9580         0.6850    98.5      62.7    375
power<-powerSFT<- 19.5

#########################################################################
#Build WGCNA net from blockwiseModules() function for GlobalNetworkPlots
minModSize=20 #for protein; in RNA, this is still 25
enforceMMS=TRUE

net <- blockwiseModules(as.matrix(t(cleanDat)),power=power,deepSplit=4,minModuleSize=minModSize,
                        mergeCutHeight=0.07, TOMDenom="mean", #detectCutHeight=0.9999,
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,reassignThresh=0.05,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=40000,replaceMissingAdjacencies = TRUE)
###I suspect that there are pairs of genes in your data that have missing values arranged such that after removing the entries that are missing in the other gene, one of the genes becomes constant.  You can get around this problem using argument
#replaceMissingAdjacencies = TRUE
#for blockwiseModules. This will instruct the underlying code to relplace missing adjacencies with zeros.
#as.matrix() to force numerioc

# Summary Module Table
table(net$colors)["grey"] #power=19 now: 110/2001

nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

# Mnum        Color Size
# turquoise      M1    turquoise  393
# blue           M2         blue  228
# brown          M3        brown  177
# yellow         M4       yellow  148
# green          M5        green  138
# red            M6          red   81
# black          M7        black   74
# pink           M8         pink   66
# magenta        M9      magenta   56
# purple        M10       purple   52
# greenyellow   M11  greenyellow   51
# tan           M12          tan   46
# salmon        M13       salmon   44
# cyan          M14         cyan   44
# midnightblue  M15 midnightblue   41
# lightcyan     M16    lightcyan   41
# grey60        M17       grey60   41
# lightgreen    M18   lightgreen   40
# lightyellow   M19  lightyellow   36
# royalblue     M20    royalblue   33
# darkred       M21      darkred   31
# darkgreen     M22    darkgreen   30

#<SKIP> (M21 is 73 members, last module--will grow) #!# My modules are all tiny, I'm going to keep them
#minModSize=73
# If necessary, return module members of small modules below size minSize=X to grey
#if (enforceMMS) {
#  removedModules<-orderedModules[which(modules<minModSize),"Color"]
#  for(i in removedModules) { net$colors[net$colors==i] <- "grey" }
#  for(i in removedModules) { net$MEs[,paste0("ME",i)] <- NULL }
#  
#  nModules<-length(table(net$colors))-1
#  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
#  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
#  modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
#  as.data.frame(cbind(orderedModules,Size=modules))
#}
#minModSize=20
#<END SKIP>

#calculate kME table up front, in case we need to correct color assignments
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
net$MEs <- MEs
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


table(net$colors)["grey"]  #still 110 grey non-module members



##ITERATIVE until condition met that all module members are at least 0.28 kMEintramodule.
#Go back and do final algorithm fix of module colors (remove kMEintramodule<0.28 members, reassign grey with kMEintramodule>0.35; max difference from kMEmax<0.10)

retry=TRUE;
kMEmaxDiff=0.1
reassignIfGT=0.30
greyIfLT=0.30
iter=1;
while (retry) {
  cat(paste0("\nkME table Cleanup, processing iteration ",iter,"..."))
  colorVecFixed<-colorVecBackup<-net$colors
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEintramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]) )  #all sig digits (no rounding), so max will be unique.
  colorVecFixed[kMEintramoduleVector<greyIfLT]<-"grey"
  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEmaxColorsVec<-apply( as.data.frame(cbind(kMEmaxVec,kMEdat)),1, function(x) gsub("kME","",colnames(kMEdat)[which(x==x[1])[2]-1]) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffTooBig<-(kMEmaxVec-kMEintramoduleVector) >= kMEmaxDiff
  colorVecFixed[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )] <- kMEmaxColorsVec[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )]
  net$colors<-colorVecFixed
  
  #  table(net$colors)["grey"]  #decreased to 2,771/8604
  #decreased to 2,771/8,587;   ROSMAP 382: 2,350/8816
  
  # Are colors still in rank order? -- put them in order by recoloring modules that changed rank
  sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]
  
  oldcolors <- names(sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"])
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==oldcolors[i]]<-paste0("proxy",labels2colors(i))
  }
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==paste0("proxy",labels2colors(i))]<-labels2colors(i)
  }
  
  # one can check that colors are in order by size now
  #sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]
  
  # recalculate kME table, since we have corrected color assignments
  MEs<-tmpMEs<-data.frame()
  MEList = moduleEigengenes(t(cleanDat), colors = net$colors, verbose=0)
  MEs = orderMEs(MEList$eigengenes)
  net$MEs <- MEs
  colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
  rownames(MEs)<-rownames(numericMeta)
  
  tmpMEs <- MEs #net$MEs
  colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
  MEs[,"grey"] <- NULL
  tmpMEs[,"MEgrey"] <- NULL
  
  kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")
  
  # recheck min kMEintramodule and max diff from kMEmax
  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEsIntramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { 1 } ) #grey proteins set to dummy value of 1 (ignore)
  
  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffCalc<- kMEmaxVec-kMEintramoduleVector
  if (min(kMEsIntramoduleVector)>=greyIfLT & max(kMEmaxDiffCalc)<=kMEmaxDiff) { cat(paste0("\nkME table 'clean' in ",iter," iterations.")); retry=FALSE; }
  iter=iter+1
  if (iter>30) break; #**
}
#** breaks after iteration 30 if did not reach criteria.


nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

#power=19.5; 11 iterations clean:
# turquoise      M1    turquoise  394
# blue           M2         blue  226
# brown          M3        brown  179
# yellow         M4       yellow  154
# green          M5        green  144
# red            M6          red   79
# black          M7        black   71
# pink           M8         pink   62
# magenta        M9      magenta   61
# purple        M10       purple   59
# greenyellow   M11  greenyellow   57
# tan           M12          tan   55
# salmon        M13       salmon   55
# cyan          M14         cyan   49
# midnightblue  M15 midnightblue   48
# lightcyan     M16    lightcyan   47
# grey60        M17       grey60   46
# lightgreen    M18   lightgreen   44
# lightyellow   M19  lightyellow   44
# royalblue     M20    royalblue   42
# darkred       M21      darkred   39
# darkgreen     M22    darkgreen   35

#recalculate kME table, since we have corrected color assignments
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
net$MEs <- MEs
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")

FileBaseName="Fib.Sen.WGCNA.MMS20.pwr19.5.lol"

## Write Module Membership/kME table!
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
colnames(kMEtable)[1:2]<-c("Orig Order","UniqueID")
write.table(kMEtable,file=paste0(outputtabs,"/ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)

## Prepare traits for Grouping
numericMeta.original<-numericMeta

table(numericMeta$treatment)
# control   copper peroxide 
# 3        6        6

table(numericMeta$serum1)
# 0 1  0 = FALSE
# 9 6

serum1<-cuso4<-h2o2<-rep(NA,nrow(numericMeta))
serum1[which(numericMeta$treatment=="control")]<-cuso4[which(numericMeta$treatment=="control")]<-h2o2[which(numericMeta$treatment=="control")]<-0
h2o2[which(numericMeta$treatment=="peroxide")]<-1
cuso4[which(numericMeta$treatment=="copper")]<-1
serum1[which(numericMeta$serum=="1")]<-1
serum1[which(numericMeta$serum=="0")]<-0
numericMeta<-numericMeta.original
numericMeta<-data.frame(treatment=numericMeta$treatment,
                        peroxide=h2o2,
                        copper=cuso4,
                        serum1=serum1,
                        treatment.serum1=paste0(numericMeta$treatment,".",numericMeta$serum))
rownames(numericMeta)<-numericMeta.original$sample
Grouping <- numericMeta$treatment
outputSuffix=FileBaseName
save.image(paste0(".saved.image.",FileBaseName,".Rdata"))  #overwrites

###################################################################################################################################
# Global Network Plots
###################################################################################################################################

library(Cairo)
for (cairoPass in 0:1) {
  if (cairoPass==0) pdf(file=paste0(outputfigs,"/GlobalNetPlots-",outputSuffix,".pdf"),width=16,height=12)
  if (cairoPass==1) CairoPDF(file=paste0(outputfigs,"/GlobalNetPlots-",outputSuffix,"-Cairo.pdf"),width=16,height=12)
  
  ## Plot dendrogram with module colors and trait correlations
  numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
  
  geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
  rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
  geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
  rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]
  
  geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
  rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
  geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
  rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]
  
  ## Page 1 -- Protein level WGCNA Dendrogram (can be left out -- if very hard to load in Acrobat)
  plotDendroAndColors(dendro=net$dendrograms[[1]],
                      colors=t(rbind(net$colors,geneSigColors)),
                      cex.dendroLabels=1.2,addGuide=TRUE,
                      dendroLabels=FALSE,
                      groupLabels=c("Module Colors",colnames(numericMeta)[numericIndices]))
  
  ## Page 2 -- Plot eigengene dendrogram/heatmap - using bicor
  plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4.6,2,0.3), marDendro = c(0,4,2,2.4),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))
  
  ## Find differences between treatment Groups (e.g. copper, peroxide, serum)
  regvars <- data.frame(as.factor(numericMeta[,"treatment"]),as.factor(numericMeta[,"peroxide"]),as.factor(numericMeta[,"copper"]),as.factor(numericMeta[,"serum1"]))
  colnames(regvars) <- c("treatment","peroxide","copper","serum1") ## data frame with covaraites in case we want to try multivariate regression
  #aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
  lm1 <- lm(data.matrix(MEs)~treatment,data=regvars)
  
  pvec <- rep(NA,ncol(MEs))
  for (i in 1:ncol(MEs)) {
    f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
    pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
  }
  names(pvec) <- colnames(MEs)
  
  ## Find differences between DxGroup.Sex Groups
  regvars <- data.frame(as.factor(numericMeta[,"treatment.serum1"]),as.factor(numericMeta[,"serum1"]))
  colnames(regvars) <- c("treatment","serum1") ## data frame with covaraites in case we want to try multivariate regression
  #aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
  lm1 <- lm(data.matrix(MEs)~treatment,data=regvars)  
  pvec.treatment.serum1 <- rep(NA,ncol(MEs))
  for (i in 1:ncol(MEs)) {
    f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
    pvec.treatment.serum1[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
  }
  names(pvec.treatment.serum1) <- colnames(MEs)  
  
  
  ## Page 3 -- Plot eigengene-trait correlations - using bicor
  library(RColorBrewer)
  MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
  moduleTraitCor <- MEcors$bicor
  moduleTraitPvalue <- MEcors$p
  
  
  textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
  #textMatrix = paste(signif(moduleTraitCor, 2), " (",
  #  signif(moduleTraitPvalue, 1), ")", sep = "");
  #dim(textMatrix) = dim(moduleTraitCor)
  par(mfrow=c(1,1))
  par(mar = c(6, 8.5, 3, 3));
  
  ## Display the correlation values within a heatmap plot
  colvec <- rep("white",1500)
  colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
  colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
  labeledHeatmap(Matrix = apply(moduleTraitPvalue,2,as.numeric),
                 xLabels = colnames(numericMeta)[numericIndices],
                 yLabels = paste0("ME",names(MEs)),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = colvec,
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(0,0.15),
                 main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
                 cex.main=0.8)
  
  
  ######################
  ## Page 4 -- Plot eigengene-trait heatmap custom - using bicor
  
  numericMetaCustom<-numericMeta[,sort(numericIndices)]
  MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
  moduleTraitCor <- MEcors$bicor
  moduleTraitPvalue <- MEcors$p
  
  moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
  moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")
  
  textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  #textMatrix = gsub("()", "", textMatrix,fixed=TRUE)
  
  labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
  labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
  for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
  for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }
  
  #rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row
  
  par(mar=c(16, 12, 3, 3) )
  par(mfrow=c(1,1))
  
  bw<-colorRampPalette(c("#0058CC", "white"))
  wr<-colorRampPalette(c("white", "#CC3300"))
  
  colvec<-c(bw(50),wr(50))
  
  labeledHeatmap(Matrix = t(moduleTraitCor)[,],
                 yLabels = colnames(numericMetaCustom),
                 xLabels = labelMat[,2],
                 xSymbols = labelMat[,1],
                 xColorLabels=TRUE,
                 colors = colvec,
                 textMatrix = t(textMatrix)[,],
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
                 verticalSeparator.col = 1,
                 verticalSeparator.lty = 1,
                 verticalSeparator.lwd = 1,
                 verticalSeparator.ext = 0,
                 horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
                 horizontalSeparator.col = 1,
                 horizontalSeparator.lty = 1,
                 horizontalSeparator.lwd = 1,
                 horizontalSeparator.ext = 0,
                 zlim = c(-1,1),
                 main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
                 cex.main=0.8)
  
  
  ## Pages 5 & 6 -- Plot clustered heatmaps (clustering MEs, pg 5; MEs and Samples, pg 6) -- annotate important metadata, plot the eigengenes!
  toplot <- MEs
  Grouping <- numericMeta$treatment
  serum1=numericMeta$serum1
  serum1[serum1==0]<-"no serum"
  serum1[serum1==1]<-"serum"
  serum1<-factor(serum1)
  
  metdat <- data.frame(Status=Grouping,serum1=serum1)
  
  colnames(toplot) <- colnames(MEs)
  rownames(toplot) <- rownames(MEs)
  toplot <- t(toplot) #may have to transpose sometimes: t(toplot)
  
  rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  | 1-way AOV P.overall=",signif(pvec,2),sep="")
  
  ## plot heatmaps of eigenproteins
  heatmapLegendColors=list('Status'=c("yellow2","darkturquoise","red","darkorange"), ## ADAD    CT  LOAD sEOAD    (alphabetical)
                           'serum1'=c("white","seagreen3"),
                           'Modules'=sort(colnames(MEs)))
  
  
  toplot.winsorized<-toplot
  minRangeExtent=min(abs(range(toplot)))
  toplot.winsorized[toplot < -minRangeExtent] <- -minRangeExtent
  toplot.winsorized[toplot >  minRangeExtent] <-  minRangeExtent
  
  
  library(NMF)
  par(mfrow=c(1,1))
  aheatmap(x=na.omit(toplot.winsorized), ## Numeric Matrix
           main="Plot of Eigengene-Trait Relationships - SAMPLES IN BATCH ORDER (MEs Winsorized)",
           annCol=metdat,
           annRow=data.frame(Modules=colnames(t(na.omit(t(MEs))))),
           annColors=heatmapLegendColors,
           border=list(matrix = TRUE),
           scale="none",
           distfun="correlation",hclustfun="average", ## Clustering options
           cexRow=0.8, ## Character sizes
           cexCol=0.8,
           col=blueWhiteRed(100), ## Color map scheme
           treeheight=80,
           Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order
  aheatmap(x=na.omit(toplot.winsorized), ## Numeric Matrix
           main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED (MEs Winsorized)",
           annCol=metdat,
           annRow=data.frame(Modules=colnames(t(na.omit(t(MEs))))),
           annColors=heatmapLegendColors,
           border=list(matrix = TRUE),
           scale="none",
           distfun="correlation",hclustfun="average", ## Clustering options
           cexRow=0.8, ## Character sizes
           cexCol=0.8,
           col=blueWhiteRed(100), ## Color map scheme
           treeheight=80,
           Rowv=TRUE,Colv=TRUE) ## Cluster columns
  #dev.off()
  
  
  
  ## Get module-trait bicor correlations (append to verboseScatterplot title below)
  numericMetaCustom<-numericMeta[,sort(numericIndices)]
  MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
  moduleTraitCor <- MEcors$bicor
  moduleTraitPvalue <- MEcors$p
  
  
  #CairoPDF(file=paste0(outputfigs,"/GlobalNetPlots(BoxPlots)_",FileBaseName,"-CAIRO.pdf"),width=18,height=11.25)
  ##pdf(file=paste0(outputfigs,"/GlobalNetPlots(BoxPlots)_",FileBaseName,".pdf"),width=18,height=11.25)
  
  #par(mfrow=c(4,6))
  par(mar=c(4.5,6,4.5,1.5))
  
  layout(matrix(c(1:20), nrow = 4, ncol = 5, byrow=TRUE),
         heights = c(0.9,0.9,0.9,0.9,0.9), # Heights of the four rows
         widths = c(1.25,0.9,0.9,1.15,0.9)) # Widths of the 5 columns
  #sapply(c(1:20),layout.show)
  
  library(beeswarm)
  library(gplots)
  
  for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
    grey40ifBlack=if (colnames(MEs)[i]=="black") { "grey40" } else { "black" }
    
    #first row
    titlecolor<-if(signif(pvec,2)[i] <0.05) { "red" } else { "black" }
    boxplot(toplot[i,]~factor(Grouping,c("control","peroxide","copper","serum")),border=grey40ifBlack,col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\n1-way AOV p = ",signif(pvec,2)[i]),xlab=NULL,las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
    transcol=paste0(col2hex(colnames(MEs)[i]),"99")
    beeswarm(toplot[i,]~factor(Grouping,c("control","peroxide","copper","serum")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col=grey40ifBlack,cex=0.8,corral="gutter") #more like prism
    
    this.tTest=t.test(toplot[i,which(Grouping=="control")], toplot[i,which(Grouping=="peroxide")], var.equal=TRUE, alternative="two.sided")$p.value
    titlecolor<-if(signif(this.tTest,2) <0.05) { if (mean(toplot[i,which(Grouping=="peroxide")]) > mean(toplot[i,which(Grouping=="control")]) ) { "red" } else { "blue" }} else { "black" }
    boxplot(toplot[i,]~factor(Grouping,c("control","peroxide")),border=grey40ifBlack,col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nT-test p = ",signif(this.tTest,2)),xlab=NULL,las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
    transcol=paste0(col2hex(colnames(MEs)[i]),"99")
    beeswarm(toplot[i,]~factor(Grouping,c("control","peroxide")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col=grey40ifBlack,cex=0.8,corral="gutter") #more like prism
    
    this.tTest=t.test(toplot[i,which(Grouping=="control")], toplot[i,which(Grouping=="copper")], var.equal=TRUE, alternative="two.sided")$p.value
    titlecolor<-if(signif(this.tTest,2) <0.05) { if (mean(toplot[i,which(Grouping=="copper")]) > mean(toplot[i,which(Grouping=="control")]) ) { "red" } else { "blue" }} else { "black" }
    boxplot(toplot[i,]~factor(Grouping,c("control","copper")),border=grey40ifBlack,col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nT-test p = ",signif(this.tTest,2)),xlab=NULL,las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
    transcol=paste0(col2hex(colnames(MEs)[i]),"99")
    beeswarm(toplot[i,]~factor(Grouping,c("control","copper")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col=grey40ifBlack,cex=0.8,corral="gutter") #more like prism
    
        #third row
    titlecolor<-if(signif(pvec.treatment.serum1,2)[i] <0.05) { "red" } else { "black" }
    boxplot(toplot[i,]~factor(numericMeta$treatment.serum1,c("control.0","control.1","peroxide.0","peroxide.1","copper.0","copper.1")),border=grey40ifBlack,col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\n1-way AOV p = ",signif(pvec.treatment.serum1,2)[i]),xlab="Treatment X Serum\n(serum=1)",las=1,col.main=titlecolor)  #no outliers: ,outline=FALSE)
    transcol=paste0(col2hex(colnames(MEs)[i]),"99")
    beeswarm(toplot[i,]~factor(numericMeta$treatment.serum1,c("control.0","control.1","peroxide.0","peroxide.1","copper.0","copper.1")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col=grey40ifBlack,cex=0.8,corral="gutter") #more like prism
    
    verboseScatterplot(x=numericMeta[,"serum1"],y=toplot[i,],xlab="Serum (serum=1)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=grey40ifBlack,bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"serum1"],2),", p=",signif(moduleTraitPvalue[i,"serum1"],2),"\n"),col.main=if(moduleTraitPvalue[i,"serum1"]<0.05) { "red" } else { "black" })
    
    
    
    
  }
  
  dev.off()
}  #ends for (cairoPass in 0:1)

############################################################################################
# iGRAPHs (Multiple Toggle Options, e.g. BioGRID interactome overlap) // CONNECTIVITY PLOT #
############################################################################################

#(re)Calculate MEs (no grey)
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL


#specify filename
baseName=FileBaseName #inputExprMat #PDF filename base
#check these settings
#######################################################
PPIedges=FALSE  #TRUE will be a lot slower...
vertexsize=8   #8 for regular, 16 for large balls
species="human" #current option "mouse" will convert bioGRID to mouse symbols before drawing PPI edges.
#CAIRO=TRUE      #Open/Write PDF using CairoPDF or pdf functions; now both file versions are created
########################################################

#(re)Establish M# for module colors table
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))


#load interactome [full bioGRID 3.4 (01/25/2016)]
#bioGrid <- read.table(file=paste0(datadir,"PPIs/BIOGRID-Homo_sapiens-3.4.133.SIMPLEsymbols_editpad.txt"),sep="\t",header=TRUE)
#------


if (species=="mouse") {
  library(biomaRt)
  human = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
  mouse = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
  
  ## If converting bioGrid from one species to another
  genelist.convert<-getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=unique(c(as.vector(bioGrid[,"From"]),as.vector(bioGrid[,"To"]))), mart=human, attributesL=c("mgi_symbol","external_gene_name"),martL = mouse)
  
  bioGrid.mouse.convert<-cbind(genelist.convert[match(bioGrid[,1],genelist.convert[,1]),3],genelist.convert[match(bioGrid[,2],genelist.convert[,1]),3])
  bioGrid.human<-bioGrid
  bioGrid<-bioGrid.mouse.convert
  dim(bioGrid)
  bioGrid<-as.data.frame(bioGrid[-which(is.na(bioGrid)),])
  colnames(bioGrid)<-colnames(bioGrid.human)
  dim(bioGrid)
}


if (PPIedges==FALSE) { bioGrid<-data.frame(From=c(1),To=c(1)) } #if you want to skip bolding edges for PPIs from BioGrid

#Symbol to add to iGraphs that don't include it, if checking for protein-protein interaction edges (add 4, one for each corner; can edit output later to remove unwanted ones.)
symbols2Add=c() #c("PiB","APP","SGIP1","APOE")

#ADmagmaList<-read.csv(file="c:/Users/Eric Dammer/Documents/wgcnatest/magma/input/AD_Mean_Zstat_genePValues_1234genesFINAL.csv",header=TRUE)
RNAbindingGOlist<- c()

#Get KMEs
# KME.vis=signedKME(cleanDat, kMEdat,corFnc="bicor",outputColumnName = "");
KME.vis=kMEdat[,] # -ncol(kMEdat) #minus grey column if calculating signedKME on the fly
annot=as.data.frame(do.call("rbind",strsplit(as.character(rownames(cleanDat)),"[|]")))
KME.vis$Symbol=annot$V1
geneInfo.vis=as.data.frame(cbind(KME.vis$Symbol,net$colors, KME.vis))  

geneInfo.vis=geneInfo.vis[,-ncol(geneInfo.vis)] # check if last column is Ensembl gene id
colnames(geneInfo.vis)[1]= "Symbol"
colnames(geneInfo.vis)[2]= "Module.Color"
colnames(geneInfo.vis)<-gsub("kMEME","kME",colnames(geneInfo.vis))

library(igraph);
library(RColorBrewer);
library(WGCNA);
library(Cairo);
library(piano)
library(ontologyIndex)

softPower = power
adjacency = adjacency(cleanDat, power=softPower, type="signed",corFnc="bicor")
TOM = TOMsimilarity(adjacency)


TOM.matrix = as.matrix(TOM);
#Get the top connected genes in the module
uniquemodcolors = gsub("kME","",gsub("kMEME","",colnames(kMEdat[,]))); #-ncol(kMEdat) last column if KMEdat.vis was calculated on the fly... moduleColors
#OR SELECT MODULES INSTEAD OF ALL
#uniquemodcolors = c("orange")


for (CAIRO in c(TRUE,FALSE)) {
  
  if (CAIRO) { CairoPDF(file=paste0(outputfigs,"/iGraph_Modules-",FileBaseName,"-CAIRO-LargeNodes.pdf"),width=16,height=12) } else {
    pdf(paste0(outputfigs,"/iGraph_Modules-",FileBaseName,"-nonCAIRO-LargeNodes.pdf"),height=9,width=10)
  }
  # for (i in 1:length(sigmodcolors))  {
  # mod=sigmodcolors[i];	
  # numgenesingraph = 50;
  # numconnections2keep = 1500;
  for (mod in uniquemodcolors)  {
    #mod="darkslateblue"
    numgenesingraph = 100;
    numconnections2keep = 700;
    cat('module:',mod,'\n');
    geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="NA",]
    geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="",]
    
    colind = which(colnames(geneInfo.vis)== paste("kME",mod,sep=""));
    rowind = which(geneInfo.vis[,2]==mod);
    cat(' ',length(rowind),'probes in module\n');
    submatrix = geneInfo.vis[rowind,];
    orderind = order(submatrix[,colind],decreasing=TRUE);
    if (length(rowind) < numgenesingraph) {
      numgenesingraph = length(rowind);
      numconnections2keep = numgenesingraph/2 * (numgenesingraph/6 - 1); #added /2 and /6 9/14/2015
    }
    if (length(rowind)<10) { innercircleNum=length(rowind) } else { innercircleNum=10 }
    cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
    submatrix = submatrix[orderind[1:numgenesingraph],];
    #Identify the columns in the TOM that correspond to these hub probes
    matchind = match(submatrix$Symbol,annot$Symbol);
    reducedTOM = TOM.matrix[matchind,matchind];
    
    orderind = order(reducedTOM,decreasing=TRUE);
    connections2keep = orderind[1:numconnections2keep];
    reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] = 1;
    
    g0 <- graph.adjacency(as.matrix(reducedTOM[1:innercircleNum,1:innercircleNum]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutMata <- layout.circle(g0)
    
    if (ncol(reducedTOM) < 51 & ncol(reducedTOM) > 10) {   
      g0 <- graph.adjacency(as.matrix(reducedTOM[11:ncol(reducedTOM),11:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMatb <- layout.circle(g0)
      
      g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75)
    } else { if (ncol(reducedTOM) > 10) { #****
      
      g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMatb <- layout.circle(g0)
      
      g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMatc <- layout.circle(g0)
      g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc)
    }
      else { layoutMat <- layoutMata } #****
    }
    
    #PLOT DONE BELOW WITH ADDED PHYSICAL INTERACTION EDGES (OR USE THIS SINGLE LINE FOR NO PHYSICAL INT HIGHLIGHTING)
    #  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(rbind(submatrix$Symbol,symbol2Add)),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*16,main=paste(mod,"module"))
    
    #Add symbols2Add
    iter=as.numeric(0)
    
    if (length(symbols2Add)==4) {
      for (symbol2Add in symbols2Add) {
        iter=iter+1
        if (iter==1) { position=cbind(-0.7071068,-0.7071068) } 
        if (iter==2) { position=rbind(position,cbind(0.7071068,-0.7071068)) }
        if (iter==3) { position=rbind(position,cbind(-0.7071068,0.7071068)) } 
        if (iter==4) { position=rbind(position,cbind(0.7071068,0.7071068)) }
        #ADD APP to graph g3 (to be used if APP not already in graph)
        g2 <- vertex(1) #,size=40,color="green", label=symbol2Add
        #layoutMat <- rbind(layoutMata*0.25,rbind(layoutMatb*0.75,cbind(-0.7071068,-0.7071068)))
        g3 <- g1+g2
        V(g3)$name <- c(1:(numgenesingraph+1))
        #WORKS# plot(g3,edge.color="grey",vertex.color=c(rep(mod,nrow(layoutMat)-1),"green"),vertex.label=as.character(c(submatrix$Symbol,symbol2Add)),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=c(submatrix[,colind]^2*16,15),main=paste(mod,"module"))
        
        g1<-g3
        numgenesingraph=numgenesingraph+1
      }
    }
    
    #x***moved out of loop
    if (length(symbols2Add)==4) {  if (ncol(reducedTOM) < 51) { layoutMat <- rbind(layoutMata*0.25,rbind(layoutMatb*0.75,position)) } else { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, rbind(layoutMatc,position*1.33)) }
    } else { if (ncol(reducedTOM) < 51 & ncol(reducedTOM) > 10) { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75) } else { if(ncol(reducedTOM) > 10) { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc) } } #do not add 4 node positions if symbols2Add is empty (length 0)
    }
    symbolList<-c(as.character(submatrix$Symbol),symbols2Add)
    if (length(symbols2Add)==4) { vertexColorVec<-c(rep(mod,numgenesingraph)-4,"green","darkgreen","steelblue","darkslateblue"); vertexSizeMat<-c(submatrix[,colind]^2*vertexsize,rep(15,4)); } else { vertexColorVec<-rep(mod,numgenesingraph); vertexSizeMat<-submatrix[,colind]^2*vertexsize; }
    
    
    
    #FIND EDGES OVERLAPPING WITH BIOGRID PAIRS
    listboldEdges<-matrix(ncol=2,nrow=0)
    if(nrow(bioGrid)>0) {
      for(i in 1:numgenesingraph) {
        for(j in i:numgenesingraph) {
          if(!(length(which(bioGrid$From==symbolList[i]&bioGrid$To==symbolList[j]))+length(which(bioGrid$From==symbolList[j]&bioGrid$To==symbolList[i])))==0) { listboldEdges <- rbind(listboldEdges,c(i,j)) }
        } }
      
      ##Remove self-loops
      if (PPIedges) { if (length(nrow(listboldEdges))>0 & nrow(listboldEdges)>0) { 
        for(i in nrow(listboldEdges):1) {
          if (i<=nrow(listboldEdges)) {
            if(listboldEdges[i,1]==listboldEdges[i,2]) {
              listboldEdges <-listboldEdges[-i,]
              if (!length(nrow(listboldEdges))>0) { 
                cat('NO PHYSICAL INTERACTIONS FOUND FOR THIS MODULE.\n')
                break #HANDLE NO BOLD EDGES CASE
              }
            }
          }
        }
      }
      }
      
    }
    
    if (is.vector(listboldEdges)) { listboldEdges<-t(data.frame(listboldEdges)) }
    
    newgraph <- g1 %>%
      #  delete_edges(listboldEdges) %>%
      set_edge_attr("color",value="lightgrey") %>%
      set_edge_attr("width",value=1) %>%
      set_edge_attr("curved",value=0) 
    
    
    if (length(symbols2Add)==4) {
      if (is.vector(listboldEdges)) { 
        cat('ONLY 1 SINGLE INTERACTION FOUND.\n')
        newgraph <- newgraph %>%
          add_edges(c(listboldEdges[1],listboldEdges[2]),color="steelblue",width=2,curved=0)
      } else {
        if (dim(listboldEdges)[1]>0) { 
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[1]))==0)&(listboldEdges[k,1]==(numgenesingraph-3)|listboldEdges[k,2]==(numgenesingraph-3))) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="#33BB33",width=3,curved=0)
            }
          }
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[2]))==0)&(listboldEdges[k,1]==(numgenesingraph-2)|listboldEdges[k,2]==(numgenesingraph-2))) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="#338833",width=3,curved=0)
            }
          }
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[3]))==0)&(listboldEdges[k,1]==(numgenesingraph-1)|listboldEdges[k,2]==(numgenesingraph-1))) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="steelblue",width=3,curved=0)
            }
          }
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[4]))==0)&(listboldEdges[k,1]==numgenesingraph|listboldEdges[k,2]==numgenesingraph)) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="darkslateblue",width=3,curved=0)
            }
          }
          
        }
      }
    } else {
      if (dim(listboldEdges)[1]>0) { 
        for(k in 1:nrow(listboldEdges)) {
          newgraph <- newgraph %>%
            add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="#33BB33",width=3,curved=0)
        }
      }
    }
    
    highlightNodes<-match(RNAbindingGOlist,symbolList) #ADmagmaList$Gene.Symbol,symbolList
    highlightNodes<-sort(na.omit(highlightNodes))
    if(length(highlightNodes)>0) { if(mod=="yellow") {vertexColorVec[highlightNodes]<-"cyan" } else {vertexColorVec[highlightNodes]<-"yellow" } }
    
    plot(newgraph,vertex.color=vertexColorVec,vertex.label=as.character(symbolList),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout=layoutMat,vertex.size=vertexSizeMat,main=paste0(orderedModules[which(orderedModules[,"Color"]==mod),"Mnum"]," ",mod," module"))
  }
  dev.off();
  
} # for (CAIRO in c(TRUE,FALSE)) loop... <repeat>

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
Grouping=numericMeta$treatment                 # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.
parallelThreads=3                              # number of CPU threads to speed calculation (recommended, 2 or more):
NETcolors=net$colors                           # list net with slot/vector containing module color assignments; length of vector must be equal to number of rows in cleanDat.
twoGroupCorrMethod="BH"                        # default method for full FDR correction when only 2 groups present is Benjamini-Hochberg; see p.adjust(..., methods= ) options.
outputCSV=TRUE                                 # Output Full Table of statistics?  TRUE/FALSE
outFilePrefix=""                              # typically "4", or step # in pipeline being run; for output file sorting by filename.
outFileSuffix="fib.senescence"                 # A description of the project, used as a filename suffix
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
selectComps=c(3,4)          # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
flip=c(3,4)              # ANOVAout column index numbers for p values in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
# As a general rule, the group with less severe effects is usually set to be the denominator (represented by what is right of '-' in ANOVAout column names)
signifP=0.05               # p value threshold for counting Differential Expression points
useNETcolors=TRUE          # use module colors saved to ANOVAout, if available; otherwise when FALSE, specify downColor upColor, and NCcolor (must be valid R color specifications in quotes)
downColor="royalblue"      # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
upColor="red"              # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
NCcolor="grey"             # points not significant are this color if useNETcolors=FALSE
splitColors=TRUE           # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
highlightGeneProducts=c()  # c("APP|P05067","MAPT|P10636","APOE|P02649") ; a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
labelHighlighted=FALSE     # if true, highlighted spots get text labels with their rownames from ANOVAout
symbolsOnly=FALSE          # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
labelTop=0                 # maximum p below which to label all points in the PDF output; OR an integer number of top ranked most significant points to label
labelSize=4.5              # text label font size, if any labels are found (when labelHighlighted=TRUE or labelTop>0)
sameScale=FALSE            # When multiple plots are drawn, should they all be the same scale with min and max x and y ranges?
HTMLout=TRUE               # output interactive HTML copies that can be opened in browser. Requires plotly package.
outFilePrefix=""          # typically the step # in the pipeline being run
outFileSuffix="fib.senescence"
# A description of the project, used as a filename suffix
outputfigs=outputfigs        # Location to save figure file output(s)

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
options(stringsAsFactors=FALSE)

## Sample pipeline outputs required to run as inputs when modulesInMemory=TRUE and ANOVAgroups=FALSE below:
#  -- contains empty cleanDat (only ordered rownames matching net$colors are needed here)
#load("c:/BaderLabGO/PipelineSample-MEGATMT-WGCNA-NatNeurosci2022.RData")

## Sample pipeline outputs required to run as inputs when ANOVAgroups=TRUE below:
#load("c:/BaderLabGO/PipelineSample-MEGATMT-ANOVA+Volcano-NatNeurosci2022.RData")


######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#inputFile <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#inputFile <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
#INPUT CSV FILE - in the filePath folder.
#Can be formatted as Kme table from WGCNA pipeline, or
#can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#in this case, the longest list is used as background or the "universe" for the FET contingencies
#  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- ouputfigs   #gsub("//","/",outputfigs)
#Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.

outFilename <- "fib.senescence.GO"
#SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.

outputGOeliteInputs=FALSE
#If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.

maxBarsPerOntology=5
#Ontologies per ontology type, used for generating the PDF report; does not limit tabled output

GMTdatabaseFile="~/R/Senescence/Human_GO_AllPathways_with_GO_iea_December_01_2022_symbol.gmt"   # e.g. "Human_GO_AllPathways_with_GO_iea_June_01_2022_symbol.gmt"
# Current month release will be downloaded if file does not exist.
# **Specify a nonexistent file to always download the current database to this folder.**
# Database .GMT file will be saved to the specified folder with its original date-specific name.
#path/to/filename of ontology database for the appropriate species (no conversion is performed)
#BaderLab website links to their current monthly update of ontologies for Human, Mouse, and Rat, minimally
#http://download.baderlab.org/EM_Genesets/current_release/
#For more information, see documentation:  http://baderlab.org/GeneSets

panelDimensions=c(3,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches

color=c("darkseagreen3","lightsteelblue1","lightpink4","goldenrod","darkorange","gold")
#colors respectively for ontology Types:
#"Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSig.C2"
#must be valid R colors

modulesInMemory=TRUE
#uses cleanDat, net, and kMEdat from WGCNA systems biology pipeline already run, and these variables must be in memory
#inputFile will be ignored
ANOVAgroups=FALSE
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#inputFile will be ignored

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################

parallelThreads=3

removeRedundantGOterms=TRUE
#if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
GO.OBOfile<-"~/R/Senescence/go.obo"
#only used and needed if above flag to remove redundant GO terms is TRUE.
#Download from http://current.geneontology.org/ontology/go.obo will commence into the specified folder if the specified filename does not exist.
#Does not appear to be species specific, stores all GO term relations and is periodically updated.

cocluster=TRUE
#If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

######################## END OF PARAMETER VARIABLES ###################################################################################

source("GOparallel-FET.R")
GOparallel()  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.
# priority is given to modulesInMemory


#Error in the final heatmap
#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [2] not equal to array extent
#In addition: Warning messages:
#  1: closing unused connection 17 (<-NEURO-2297-PC.Eu.Emory.Edu:11273) 
#2: closing unused connection 16 (<-NEURO-2297-PC.Eu.Emory.Edu:11273) 
#3: closing unused connection 15 (<-NEURO-2297-PC.Eu.Emory.Edu:11273) 
#4: closing unused connection 14 (<-NEURO-2297-PC.Eu.Emory.Edu:11273) 
#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [2] not equal to array extent
#3.
#`colnames<-`(`*tmp*`, value = generate_dimnames(labCol, ncol(mat), 
#                                                colnames(mat)))
#2.
#aheatmap(x = data, main = "Co-clustering with manhattan distance function, ward metric", 
#         annRow = myRowAnnotation, annColors = heatmapLegendColors, 
#         border = list(matrix = TRUE), scale = "none", distfun = "manhattan", 
#         hclustfun = "ward", cexRow = 0.8, cexCol = 0.8, col = colvec, ... at GOparallel-FET.R#793
#         1.
#         GOparallel()
         



