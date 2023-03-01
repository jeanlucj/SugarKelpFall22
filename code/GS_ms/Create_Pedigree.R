# Re do the pedigree matrix

# save(GPsequenced,kelpNameColumns,SPGPs,file="What_has_been_Used_from_Step1_Create_Matrix.Rdata")

setwd("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/ReorderPedigree")
library(rrBLUP)
load("what_has_been_used_from_step1_create_matrix.Rdata")
  # A."GPsequenced": the list of GPs sequenced in 3 plates, remained 270 of them
  # B."kelpNameColumns": list of femaPar x malePar, GP in making the SP crosses had photo score >1. 244 plots
  # C.SPGPs: the list of GPs came from 2019-SP plot level, produced biomass in 2020

GP_Add<-read.csv("AddGP_to_Pedi.csv",sep=",",header=T)
  # D. "GP_Add":The additional list of GPs to add (had biomass in UCONN stock, in addition to A.B.C)

colnames(GP_Add)<-"WBiomass_Name_InPheno"
  head(GP_Add)
  dim(GP_Add)

load("FarmCPU_GAPIT.Rdata")
  # E. "geno": the marker matrix for genotyped 2018 fnders
geno2<-geno[,-1]
rownames(geno2)<-geno$taxa
  dim(geno2)  #125  
rownames(geno2) <- paste0(substring(rownames(geno2), first=1, last=2), "18", substring(rownames(geno2), first=3))
  ls()
  
source("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/is.square.matrix.R")
source("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/is.positive.definite.R")
source("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/is.symmetric.matrix.R")


###1 GP genotyped,crossed and had SP phenotypic data in dataNHpi
###2 GP not genotyped?,crossed and had SP phenotypic data in dataNHpi
###3 GP genotyped,not crossed and no SP phenotypic data in dataNHpi
###4 GP from 2019S_plot level
###5 other GP had biomass,not genotyped and not crossed

  str(GPsequenced)  # 270
  str(kelpNameColumns)  # 244
  str(unique(GPsequenced)) #270
  str(unique(kelpNameColumns$femaPar))  #106
  str(unique(kelpNameColumns$malePar))  #121

##1 
FG_seq<-unique(kelpNameColumns$femaPar[kelpNameColumns$femaPar%in%GPsequenced]) #68

MG_seq<-unique(kelpNameColumns$malePar[kelpNameColumns$malePar%in%GPsequenced]) #91

GP_seq<-c(FG_seq,MG_seq)   # a Unique list for GPs in 1
##2
FG_Noseq<-unique(kelpNameColumns$femaPar[!kelpNameColumns$femaPar%in%GPsequenced]) #38

MG_Noseq<-unique(kelpNameColumns$malePar[!kelpNameColumns$malePar%in%GPsequenced]) #30

GP_Noseq<-c(FG_Noseq,MG_Noseq)  # a Unique list for GPs in 2

##3  
GP_noPheno<-GPsequenced[!GPsequenced%in%c(kelpNameColumns$femaPar,kelpNameColumns$malePar)] #111
           # a Unique list for GPs in 3

##4 
S_GP<-SPGPs$GametophyteID

# I. find all 1,2,3,4 their fnders, order, unique the list
# II. GPs themselves, 1,2 and 3
# III.The pair of SP plots in kelpNamesColumns femaParx malePar
# IV. The rest of the SPGPs$SPcrosses (parents for 4); adding that plot for SPGPs$SPCrosses that is not in the kelpNamesColumns femaPar x malePar
# V. S_GPs for 4 themselves

find_fnder_for_GP<-function(GPstring){
  fndr<-strsplit(GPstring,split="-", fixed=T) 
  fnder_unique<-unique(sapply(fndr,function(vec) paste(vec[1:3], collapse="-")))
  return(fnder_unique)
}

#1
GP_seq_fndr<-find_fnder_for_GP(GP_seq) #64
#2
GP_Noseq_fndr<-find_fnder_for_GP(GP_Noseq) # 44
  str(GP_seq_fndr)  # 1
  str(GP_Noseq_fndr) # 2
  str(unique(c(GP_seq_fndr,GP_Noseq_fndr)))  # 70 
#3
GP_noPheno_fndr<-find_fnder_for_GP(GP_noPheno)  #61 
  str(GP_noPheno_fndr) # 3
  str(unique(c(unique(c(GP_seq_fndr,GP_Noseq_fndr)),GP_noPheno_fndr)))  #93

library(stringr)
find_fnder_for_SGP<-function(S_GPCross){
  F_M_GPs<-str_split_fixed(string=as.character(S_GPCross), "x", 2)
  S_FGP<-F_M_GPs[,1]
  S_MGP<-F_M_GPs[,2]
  S_FGP_fndr<-find_fnder_for_GP(S_FGP)  #11
  S_MGP_fndr<-find_fnder_for_GP(S_MGP)  #9
  
  return(list(S_FGP=S_FGP,S_MGP=S_MGP,S_FGP_fndr=S_FGP_fndr,S_MGP_fndr=S_MGP_fndr))
}


S_GP_its_fnders<-c(find_fnder_for_SGP(SPGPs$SPCross)$S_FGP_fndr,find_fnder_for_SGP(SPGPs$SPCross)$S_MGP_fndr)
  str(unique(c(unique(c(unique(c(GP_seq_fndr,GP_Noseq_fndr)),GP_noPheno_fndr)),S_GP_its_fnders))) # 93 

S_GP_its_parentalGP<-unique(c(find_fnder_for_SGP(SPGPs$SPCross)$S_FGP,find_fnder_for_SGP(SPGPs$SPCross)$S_MGP))


###!!! 5. Adding GPs has biomass, but not in the field data (photo score >1) nor in the genotyped list
GP_Add_unique<-unique(as.character(GP_Add$WBiomass_Name_InPheno)) # 105

# find their founders
# themselves
BiomassGP_fndr<-find_fnder_for_GP(GP_Add_unique)  

#### List of fnders for 1,2,3,4
All_fnders1<-c(GP_seq_fndr,GP_Noseq_fndr,GP_noPheno_fndr,S_GP_its_fnders)  # 189
All_fnders1_unique<-unique(All_fnders1)  # 93
  str(BiomassGP_fndr[!BiomassGP_fndr%in%All_fnders1_unique])

  #For 4, all S_GP_its_fnders were already included
  #For 1, to row64, For 2 row65 to row70, For 3 row71 to row93)
  #For 5, row 94 to row 104
  
All_fnders2<-c(All_fnders1,BiomassGP_fndr) #249

#### List of fnders for 1,2,3,4,5
All_fnders_unique<-unique(All_fnders2)   #104
  str(All_fnders_unique)                    

## I. Ped for fnders

## Use Input file E. geno2to Re-order fnders
geno2<-geno2[rownames(geno2)%in%All_fnders_unique,]  # 58
  #geno2<-geno2[rownames(geno2)%in%CrossedSP,] This only gave 47 fndrs being "both genotyped, and crossed to make SP <PhotoScore>1)
  # # The fndrs list previously used, only kept 47 of them that had genotypic data +photo score >1. In "CovList_3_As_0112_2021.Rdata"
  ##### This will be ready for outCovComb
fndrMrkData<-geno2
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T,shrink=TRUE) ## Add shrink per Deniz
fndrMrkDataImp <- mrkRelMat$imputed
mrkRelMat <- mrkRelMat$A
fndrsA<-round(mrkRelMat,digits=5)  # off diagonal decimal points not the same -> not symmetric matrix   
  diag(fndrsA) <- diag(fndrsA) + 1e-5
is.positive.definite(fndrsA)

#### This is the new order
  # a) fndrs genotyped
  # b) for the list of GPs genotyped, find their fndrs
  # 1) common fndrs in a) and b)
  # 2) fndrs in a)-1)
  # 3) fndrs in b) -a)
  # 4) the rest of the other fndrs
  
  a<-rownames(fndrsA)
  b<-find_fnder_for_GP(GPsequenced)
  
  a_b_common<-intersect(a,b) # 1. fndr+GP genotyped
  a_1<-setdiff(a,a_b_common) # 2. fndr genotyped, its GP not genotyped
  b_a<-setdiff(b,a)    # 3. its GP genotyped, fndr not genotyped
  
  nfdr<-length(All_fnders_unique)   # length of fnders !!!
  other<-setdiff(All_fnders_unique,c(a_b_common,a_1,b_a))
  
  #1:56, 57:59,60:93,94:104
  fndrs_NewOrder<-c(a_b_common,a_1,b_a,other)
  sum(fndrs_NewOrder%in%All_fnders_unique)
  
  # biphasic_NewOrder<-c(fndrs_NewOrder,All_fnders_unique[(nfdr+1):nrow(biphasicPedNH)])
  # identical(biphasic_NewOrder[(nfdr+1):length(biphasic_NewOrder)],All_fnders_unique[(nfdr+1):nrow(biphasicPedNH)])
  # 
  # biphasicPedNH_fndrs_NewOrder<-biphasicPedNH[match(biphasic_NewOrder,All_fnders_unique),]
  
  All_fnders_unique<-All_fnders_unique[order(match(All_fnders_unique,fndrs_NewOrder))]
  identical(All_fnders_unique,fndrs_NewOrder)

## Create the Ped format  
nFounders<-length(All_fnders_unique) 
fndRow <- 1:nFounders
names(fndRow) <- All_fnders_unique
fndPed <- cbind(fndRow, 0, 0)   

All_GPs123<-c(GP_seq,GP_Noseq,GP_noPheno)   # 338

## GP_Add_unique: 5. the Added_GPs that's not in genotyped nor phenotyped but had biomass  
GP_Add_unique<-GP_Add_unique[!GP_Add_unique%in%All_GPs123]  # 105->101

All_GPs<-c(All_GPs123,GP_Add_unique)   # 439
  str(All_GPs)      # 439
  length(GP_seq)     # row105 to row263 length(c(All_fnders_unique,GP_seq))
  length(GP_Noseq)   # row264 to row331 /length(c(All_fnders_unique,GP_seq,GP_Noseq))
  length(GP_noPheno)  #row332 to row442 /length(c(All_fnders_unique,GP_seq,GP_Noseq,GP_noPheno))
  length(GP_Add_unique)   #row443 to row 543 /length(c(All_fnders_unique,GP_seq,GP_Noseq,GP_noPheno,GP_Add_unique))

sum(S_GP_its_parentalGP %in%(All_GPs))==length(S_GP_its_parentalGP)   #All the GPs used for 4 are in the All_GPs list

## II. Ped for GPs #### DO I NEED TO ORDER THEM BY FEMALE firt and then BY MALE ???
GPsRow <- nFounders + 1:length(All_GPs)
names(GPsRow) <- All_GPs
founderPar<-sapply(strsplit(All_GPs, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
GPsPed <- cbind(GPsRow, fndRow[founderPar], NA)

## III. Ped for the SP plots in 1/2, aka photo score >1 plots, aka in kelpNamesColumns femaParx malePar          
sporProg1 <- 1:nrow(kelpNameColumns) + nrow(fndPed) + nrow(GPsPed)
names(sporProg1)<-paste0(kelpNameColumns[,1],"x",kelpNameColumns[,2])
progPed1 <- cbind(sporProg1, GPsRow[kelpNameColumns[,1]], GPsRow[kelpNameColumns[,2]]) # row 432 to row 675 (431+244)           

## IV. Ped for the SP_GP's plots that were not in 1/2        
# The SPplots not listed in the kelpNameColumns combination (because their plot photo score may not pass 1)
sproProg2Name<-unique(droplevels(SPGPs$SPCross[!SPGPs$SPCross%in%names(sporProg1)]))        
sproProg2<-1:length(sproProg2Name)+nrow(fndPed)+nrow(GPsPed)+nrow(progPed1)

names(sproProg2)<-sproProg2Name
progPed2<-cbind(sproProg2,GPsRow[str_split_fixed(string=as.character(names(sproProg2)), "x", 2)[,1]], GPsRow[str_split_fixed(string=as.character(names(sproProg2)), "x", 2)[,2]])

progPed<-rbind(progPed1,progPed2)

## V. Ped for the SP_GPs themselves
spRows<-progPed[,1]
nrow(SPGPs)==length(unique(SPGPs$GametophyteID))
SP_GP_row<-1:nrow(SPGPs)+nrow(fndPed)+nrow(GPsPed)+nrow(progPed)
names(SP_GP_row)<-SPGPs$GametophyteID
SP_GP_Ped<-cbind(SP_GP_row,spRows[as.character(SPGPs$SPCross)],NA)
  
Ped_in_Order<-rbind(fndPed,GPsPed,progPed,SP_GP_Ped)  
write.csv(Ped_in_Order,"Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv")

### calculate the aMat (pedigree based relationship matrix)_ Diploid level
  
source("calcCCmatrixBiphasic.R")
#biphasicPedNH<-read.csv("Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv.csv",sep=",",header=TRUE,row.names=1)
biphasicPedNH<-Ped_in_Order

biphasicCCmat <- calcCCmatrixBiphasic(biphasicPedNH)  #### !!!!!!! Update into calcCCmatrixHaploid() ?????
rownames(biphasicCCmat) <- colnames(biphasicCCmat) <- rownames(biphasicPedNH)
aMat <- 2 * biphasicCCmat

save(fndrsA,biphasicPedNH,aMat,file="fndrsA_biphasicPedNH_fnderOrdered.Rdata")


####### Calculating the outCovComb
rm(list=ls())
load("fndrsA_biphasicPedNH_fnderOrdered.Rdata")  
  # "fndrsA": fndrs A matrix
  # "biphasicPedNH": pedigree, Newly ordered fndrs 0116_2021, GPs also already ordered by 0114_2021
  # "aMat": A matrix estimated from biphasicPedNH
load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/SugarKelpBreeding_NoGitPush/GenotypicData_for_SugarKelpBreeding/GPsAmat_NA0.8_P1P2P3_09282020.Rdata")
  # Load the GPs A matrix, estimated from A.mat()
  # Further reduced the "SL18-LD-13-Mg-3"
GPsA<-GPsA[!rownames(GPsA)=="SL18-LD-13-Mg-3",!colnames(GPsA)=="SL18-LD-13-Mg-3"]
GPsA<-round(GPsA,digits=5)   ### 2.GPs A
  diag(GPsA) <- diag(GPsA) + 1e-5
  is.positive.definite(GPsA)  # GPsA2 calculated by hand (above)  # 278

diag(aMat) <- diag(aMat) + 1e-5
is.positive.definite(aMat)     ### 3. pedigree

  identical(rownames(aMat),rownames(Ped_in_Order))
  sum(rownames(fndrsA)%in%rownames(aMat))   # 58
  
  sum(rownames(fndrsA)==colnames(fndrsA))
  sum(rownames(GPsA)==colnames(GPsA))  # 278
  sum(rownames(aMat)==colnames(aMat))
  sum(rownames(fndrsA)%in%rownames(aMat))
  sum(rownames(GPsA)%in%rownames(aMat))  # 270  !///!! This became a problem
write.csv(rownames(GPsA)[!rownames(GPsA)%in%rownames(aMat)],"GPsA_not_in_aMat.csv")

GPsA<-GPsA[rownames(GPsA)%in%rownames(aMat),colnames(GPsA)%in%rownames(aMat)]  # 270 x 270

save(fndrsA,GPsA,aMat,file="CovList_3_As_0116_2021.Rdata")

##{} Run this in terminal
#This was where old files were calculated
#setwd("/local/workdir/mh865/SNPCalling/Saccharina_latissima/bamAll/bam/mpileup/Filtering/CovComb")

setwd("/local/workdir/mh865/outCovComb/")
load("CovList_3_As_0116_2021.Rdata")

library(CovCombR)

### 1. NO initial, 3-list
CovList<-NULL
CovList[[1]]<-fndrsA ## fndrsA
CovList[[2]]<-GPsA  ## Further RMed the "SL18-LD-13-Mg-3"
CovList[[3]]<-aMat   

### 4. amat initial, 3-list, add weight
weights<-c(2,2,1)
outCovComb4<-CovComb(CovList,nu=1500,w=weights,Kinit=aMat) 

outCovComb4_dipOrder<-outCovComb4[match(rownames(aMat),rownames(outCovComb4)),match(colnames(aMat),colnames(outCovComb4))]
  # Reorder the outCovComb4 to aMat individuals!!
save(outCovComb4_dipOrder,file="outCovComb_dip_0116_2021.Rdata")
write.csv(rownames(outCovComb4),"outCovComb4_from_ordered_pedi_rownames.csv")


#### Calculating outCovComb on a haploid base, go to Step 1)haploid.R
####
####



#### FOR JL: Ignore these below part for now
#### Compare two outCovComb4_old not ordered, and the outCovComb4_new ordred. Both at diploid level, r=0.99
# In terminal
load("outCovComb_files_0112_2021.Rdata")

load("CovList_3_As_0112_2021.Rdata")
outCovComb4_order<-outCovComb4[match(rownames(aMat),rownames(outCovComb4)),match(colnames(aMat),colnames(outCovComb4))]

load("/local/workdir/mh865/SNPCalling/Saccharina_latissima/bamAll/bam/mpileup/Filtering/outCovComb4_10012020_withSGP_866.Rdata")
outCovComb4_NonOrder<-outCovComb4[match(rownames(aMat),rownames(outCovComb4)),match(colnames(aMat),colnames(outCovComb4))]

install.packages("ape")
library(ape)
mantel.test(as.matrix(outCovComb4_order),as.matrix(outCovComb4_NonOrder),trials=99)

library(ade4)
mantel.rtest(dist(as.matrix(outCovComb4_order)),dist(as.matrix(outCovComb4_NonOrder)),nrepet=99)

cor.test(c(outCovComb4_order),c(outCovComb4_NonOrder))
#cor.test(c(outCovComb4_order),c(outCovComb4_NonOrder))

# #Pearson's product-moment correlation
# 
# data:  c(outCovComb4_order) and c(outCovComb4_NonOrder)
# t = 12596, df = 749950, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9976344 0.9976557
# sample estimates:
#       cor 
# 0.9976451
# 
# mantel.test
# #$z.stat
# [1] 4287.168
# 
# $p
# [1] 0.001
# 
# $alternative
# [1] "two.sided"

#mantel.rtest
# Observation: 0.9956624
# 
# Based on 99 replicates
# Simulated p-value: 0.01
# Alternative hypothesis: greater
# 
# Std.Obs   Expectation      Variance
# 42.7492404702 -0.0030156889  0.0005457506


