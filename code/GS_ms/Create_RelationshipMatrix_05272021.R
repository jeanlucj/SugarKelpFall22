### Re-make the relationship matrix
#1. Pedigree CC matrix calculated at haploid level, condense using JL function
#2. FndrA calculated at diploid level via a.mat()
#3. GPs A calculated using the haploid function
rm(list=ls())
#write.csv(Ped_in_Order,"Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv") # This is the most updated Pedi-file

datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data"

### 1. Pedigree-> CC matrix at haploid level -> condense
# To get the pedigree at diploid level": "biphasicPedNH"
library(here)
here()
biphasicPedNH<-read.csv(here("TraitAnalyses201003/ReorderPedigree","Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv"),sep=",",header=TRUE,row.names=1)

# To get the function calculating haploid level ccMatrix, also condensing
source(here("TraitAnalyses201003/Making_haploid_CovComb","calcCCmatrixBiphasic.R"))
source(here("TraitAnalyses201003/Code_10032020","condenseMixedPloidyRelMat.R"))

# Use the diploid Pedigree file, ccMatrix change rows names
biphasicPedNH<-as.matrix(biphasicPedNH)

HaploidCCcal <- calcCCmatrixHaploid(biphasicPedNH)  ###

biphasichapccMat<-HaploidCCcal$hccMat   ###!!! 2. CCmatrix/aMat at haploid level for pedigree
hapOnePointer<-HaploidCCcal$hapOnePointer  ###!!! For condensing
dipTwoPointers<-HaploidCCcal$dipTwoPointers  ###!!! For condensing

### This is to rename the rownames out of the biphasichapccMat
# pedinames<-NULL
# for (i in 1:nrow(biphasicPedNH)){
#   if(!is.na(biphasicPedNH[i,2]) & !is.na(biphasicPedNH[i,3])){
#     # this is diploid, both parents have numbers
#     pedinames<-c(pedinames,rownames(biphasicPedNH)[i],paste0(rownames(biphasicPedNH)[i],"_2"))
#   }else{
#     # this is GP, one parent is missing
#     pedinames<-c(pedinames,rownames(biphasicPedNH)[i])
#   }
# }
#
# rownames(biphasichapccMat)<-colnames(biphasichapccMat)<-pedinames

### Condense it here: Matrix 2
Pediconden<-condenseMixedPloidyRelMat(hccMat=biphasichapccMat,hapOnePointer =hapOnePointer,dipTwoPointers=dipTwoPointers)
# 866 x 866
rownames(Pediconden)<-colnames(Pediconden)<-rownames(biphasicPedNH)

#2. Input fndrsA at diploid level files:
load(paste0(datafdr,"/CovList_3_As_0116_2021.Rdata"))  #####!!!!! Now fndrs has 58 individuals
# To get the list of fndrs from diploid level: "fndrsA"
dim(fndrsA)  # 58 x 58
fndrsA_diploid<-rownames(fndrsA)

### 3.
load(paste0(datafdr,"/Newly_saved_3_As_for_CovComb_0116_2021.Rdata"))
GPsA<-GPsA  # 269 x 269

diag(Pediconden) <- diag(Pediconden ) + 1e-5
is.positive.definite(Pediconden)
save(Pediconden,fndrsA,GPsA,file="Newly_saved_3_As_for_CovComb_0420_2021.Rdata")
#### Note GPs data  # did not redo this part on 0116_2021
                    # did not redo this part on 0420_2021

######{{}} Did this in the Terminal to estimate the RelationshipMatrix2
### 3. To get the GPs the raw SNPs data: "GPSNP" is the raw SNPs"
#load(here("../2020_2019_Phenotypic_Data/SugarKelpBreeding_NoGitPush/GenotypicData_for_SugarKelpBreeding","GPs_mainGenome_NA0.8_P1P2P3.Rdata"))
# Calculate GRM as (W %*% W^T) / sum(locusDosageVariance)
# To calculate p, the function expects dosage coding to be 0, 1, 2 for diploid
# and 0, 1 for haploid
### Calculate the relationship matrix for GPs
calcGenomicRelationshipMatrix <- function(locusMat, ploidy=2){
  if (!any(ploidy == 1:2)) stop("Ploidy must be 1 or 2")
  freq <- colMeans(locusMat) / ploidy
  locusMat <- scale(locusMat, center=T, scale=F)
  return(tcrossprod(locusMat) / sum(ploidy*freq*(1-freq)))
}

###IF to RUN in Terminal
load("/local/workdir/mh865/SNPCalling/Saccharina_latissima/bamAll/bam/mpileup/Filtering/GPs_mainGenome_NA0.8.Rdata")
#### Need to redo the formating of the raw marker data again
GPs0<-GPSNP   #
rownames(GPs0)<-paste0(GPs0$chrom,"_",GPs0$pos)
map<-GPs0[,1:2]
GPs<-GPs0[,-c(1:2)]

#### Correct the names of GPs
library(stringr)
names<-colnames(GPs)
names<-str_replace_all(names,"SA-","SA18-")
names<-str_replace_all(names,"SL-","SL18-")
names<-str_replace_all(names,"FG-","FG")
names<-str_replace_all(names,"MG-","MG")

colnames(GPs)<-names

## RM checks, the mistake one (RM SL-PI-1-FG-1,SL-CT1-FG-3)
IndiRM<-c("SL18-PI-1-FG1","3","SL18-CT1-FG3","SL18-CT1-MG2","SL18-CT1-MG3","SL18-OI-15-Female","SL18-ME-1-FG1","SL18-ME-1-MG1","SL18-ME-1-MG2")
GPs1<-GPs[,!colnames(GPs)%in%IndiRM]
GPs<-GPs1

#Impute GP marker matrix with A.mat()
tSNP<-t(GPs)
tSNP[tSNP==2]=1  # Change to haploid level dosage

#mrkMat needs to be m xn dimension
imputeToMean <- function(mrkMat){
  imputeVec <- function(mrkVec){
    meanScore <- mean(mrkVec, na.rm=T)
    mrkVec[is.na(mrkVec)] <- meanScore
    return(mrkVec)
  }
  imputedMrkMat <- apply(mrkMat, 2, imputeVec)
  return(imputedMrkMat)
}

GPsMrkImp<-imputeToMean(tSNP)
dim(GPsMrkImp)  #269 x 909749

#Change GPs marker to 0, 1 format

onerow<-GPsMrkImp
onerow[onerow> 0 & onerow< 0.5]=0
onerow[onerow>= 0.5 & onerow< 1]=1
onerow[onerow> 1]=1

### Different ways of estimating GPsAs
GPsA<-calcGenomicRelationshipMatrix(GPsMrkImp,ploidy=1)

GPsA2<-calcGenomicRelationshipMatrix(onerow,ploidy=1) ##Use This !!!## This function assumes data to be 0 and 1 for haploid

cultevo::mantel.test(dist(GPsA),dist(GPsA2),trials=99) # 0.969

M<-GPsMrkImp
GPsA3<-cov(t(M))/mean(diag(cov(t(M))))
mantel.test(dist(GPsA),dist(GPsA3),trials=99) # 0.984

M2<-onerow
GPsA4<-cov(t(M2))/mean(diag(cov(t(M))))
mantel.test(dist(GPsA2),dist(GPsA4),trials=99) # 0.9

### Choose the onerow, calcGenomicRelationshipMatrix(), so GPsA2

source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.square.matrix.R")
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.positive.definite.R")
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.symmetric.matrix.R")

GPsA<-GPsA2   #269
diag(fndrA) <- diag(fndrA) + 1e-5
is.positive.definite(fndrA)   #116
diag(GPsA) <- diag(GPsA) + 1e-5    #269
is.positive.definite(GPsA)
diag(biphasichapccMat) <- diag(biphasichapccMat) + 1e-5   #1215
is.positive.definite(biphasichapccMat)

# ### RM the extra fndr individuals that's not in the pedigree list
# RMextraFndr<-which(!rownames(fndrA)%in%rownames(biphasichapccMat))
# fndrA2<-fndrA[-RMextraFndr,-RMextraFndr]
#   dim(fndrA2)
#   fndrA2[29:30,29:30]

### check and see Is GPsA2 individuals all in the biphasichapccMat?
which(!rownames(GPsA)%in%rownames(biphasichapccMat))
GPsA<-GPsA[rownames(GPsA)%in%rownames(aMat),colnames(GPsA)%in%rownames(aMat)]  # 269 x 269
which(!rownames(fndrsA)%in%rownames(biphasichapccMat))
is.positive.definite(fndrA)

#save(fndrA,GPsA,biphasichapccMat,file="Newly_saved_3_As_for_CovComb_0116_2021.Rdata")
save(hapOnePointer,dipTwoPointers,file="Pointers_from_biphasichapccMat_0116_2021.Rdata")

###{{Done in terminal}}

#### 05_27_2021 Redo the fndrsA matrix, by multiplying to 2
rm(list=ls())
datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data"
load(paste0(datafdr,"/Newly_saved_3_As_for_CovComb_0420_2021.Rdata"))
  ls()   # Pediconde, fndrsA, GPsA
  dim(fndrsA)
  dim(GPsA)
  dim(Pediconden)
fndrsA<-fndrsA*2
save(Pediconden,GPsA,fndrsA,file=paste0(datafdr,"/Newly_saved_3_As_for_CovComb_0527_2021.Rdata"))


##{{Done in terminal}}
#cd /local/workdir/mh865/outCovComb
library(CovCombR)
load("Newly_saved_3_As_for_CovComb_0527_2021.Rdata")

### The order of fndrsA is not the first 58 of Pediconden, OK?

CovList<-NULL
CovList[[1]]<-fndrsA ## fndrsA

CovList[[2]]<-GPsA
CovList[[3]]<- Pediconden

weights<-c(2,2,1)


outCovComb4<-CovComb(Klist=CovList,nu=2000,w=weights,Kinit=Pediconden)
####### HAS TO BE RE-ORDERED !!!
outCovComb4_MixOrder<-outCovComb4[match(rownames(Pediconden),rownames(outCovComb4)),match(colnames(Pediconden),colnames(outCovComb4))]

save(outCovComb4_MixOrder,Pediconden,fndrsA,GPsA,file="outCovComb4_Mix_Conden_0527_2021.Rdata")
