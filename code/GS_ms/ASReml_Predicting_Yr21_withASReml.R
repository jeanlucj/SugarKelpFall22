#### This script contains using raw Yr19 and Yr20 data, Fitting both years' data in the model with unstructured year variance and heterogeneous error variance between years.
#### This will get the estimated h2 within each year, Genetic cor for the same trait across years, SE for those. Also the BVs of GPs being used in Yr2021
#### Yr19,Yr20 predicting Yr21 with heterogeneous error variance structure
#### Results output are in: "/local/workdir/mh865/GCA_SCA/"

rm(list=ls())
library(asreml)

#### Read in the Yr2021 data
WD2<-"/local/workdir/mh865/GCA_SCA/"
datafdr2<-paste0(WD2,"OneTime192021/data/")
load(paste0(datafdr2,"3yrs_data_plotformated_indivNot_12012021.Rdata"))    ### Included the Yr2021 SPs in pedigree
dataYr21<-droplevels(dataNHpi3yrs_C_Ash[dataNHpi3yrs_C_Ash$Year==2021,])
dataYr21GPlist<-unique(c(as.character(dataYr21$femaPar),as.character(dataYr21$malePar)))  ### Obtain the list of GPs used in making Yr21 SPs (photoscore >=2)

#### Read in Yr2019 and Yr2020 data
WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-paste0(WD,"OneTime1920/data/")
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021_866Indiv.Rdata"))   ### 866 individuals, Only included the GPs for Yr2021 in the pedigree, but not the SPs of Yr2021
load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23.rdata"))  ## Plot -- Updated Ash dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata


dataNHpi<-dataNHpiBoth_C  ##!!!!! Yr19+Yr2020 data combined
exptlSP <- as.character(dataNHpi$popChk) == "ES"
dataNHpi$group <- as.factor(ifelse(exptlSP, 1, 0))

for (col in c( "line", "block","popChk","group","Year")){ 
  dataNHpi[,col] <- factor(dataNHpi[,col])}
dataNHpi<-droplevels(dataNHpi)

### Load in Individual trait data
load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123_measurementsInfactor.rdata"))  

### Genomic relationship matrix only using Yr19+Yr20 SPs (+founders +GPs) 866*866
All_grm<-outCovComb4_MixOrder
  dim(All_grm)
  str(dataNHpi)

### Run below Plot Level VS Indi Level traits separatelly, Run only once !!! 
  
### 1. Plot Level----6 traits
  PlotTraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
  dataNH<-dataNHpi   ### !!!!! Only the plot traits
  
  traits<-PlotTraits
  for (col in traits){
    dataNH[,col]<-ifelse(dataNH[,col]==0,NA,dataNH[,col]) # Convert 0s to NAs
  }
  
### Only Run ONCE!!!!
  dataNH$wetWgtPerM0<-dataNH$wetWgtPerM
  dataNH$percDryWgt0<-dataNH$percDryWgt
  dataNH$dryWgtPerM0<-dataNH$dryWgtPerM
  dataNH$Ash0<-dataNH$Ash
  dataNH$AshFDwPM0<-dataNH$AshFDwPM
  dataNH$densityBlades0<-dataNH$densityBlades
  
  ### log transform
  for (col in traits){
    dataNH[,col]<-log(dataNH[,col]+1)
  }
  
  trtlevel<-"PlotLevel"
  traits<-PlotTraits
### Only Run ONCE!!!!
  
  
  
#### 2. Individual level to get their experimental factors.
#### The dataNHpi is already filtered for their phenotypic PHOTO SCORE (>2)
  
morTraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
  
AveIndi_to_PlotLevel<-function(dataNHpi=dataNHpi,dataNHim=dataNHim3yrs_C,popChk=TRUE){
  dataNHim<-dataNHim3yrs_C
  dataNHim$line<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "line")
  dataNHim$block<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "block")
  dataNHim$Year<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Year")
  if(popChk==TRUE){
  dataNHim$popChk<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "popChk")
            }
  dataNHim$Crosses<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Crosses")
  dataNHim$PhotoScore<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "PhotoScore")
  dataNHim<-dataNHim[which(dataNHim$PhotoScore >1),]  # 3969 rows with PhotoScore >1
   print(str(dataNHim))
   print(dim(dataNHim))
#### Note: This formating step is a MUST, IF you are seeing morTraits being factors !!!
    for (col in morTraits){
    dataNHim[,col]<-as.numeric(as.character(dataNHim[,col])) }  ##!!!! This needs to be as.character, then as.numeric !!!
####
  
####### Note: This below is especially needed if doing multi-trait genetic cor: merge Indi level to plot level
####### Averaging the individual measurements to a plot level, and merge to plot traits---> dataNH1

  library(dplyr)
  dataNHim_avg<-aggregate(cbind(bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ Year+line+block+Crosses, data=dataNHim, FUN=mean, na.action = na.omit)
  ### Plot and individual traits
  dataNHpi$MergeVar<-paste(dataNHpi$Year,dataNHpi$line,dataNHpi$block,dataNHpi$Crosses,sep="_")
  dataNHim_avg$MergeVar<-paste(dataNHim_avg$Year,dataNHim_avg$line,dataNHim_avg$block,dataNHim_avg$Crosses,sep="_")
  
  dataNH1<-merge(dataNHpi,dataNHim_avg,by.x="MergeVar",by.y="MergeVar",all.x=TRUE)
  dataNH1[is.na(as.character(dataNH1$Crosses.y)==as.character(dataNH1$Crosses.y)),]
  colnames(dataNH1)[colnames(dataNH1)=="Year.x"]<-"Year"
  colnames(dataNH1)[colnames(dataNH1)=="line.x"]<-"line"
  colnames(dataNH1)[colnames(dataNH1)=="block.x"]<-"block"
  colnames(dataNH1)[colnames(dataNH1)=="Crosses.x"]<-"Crosses"
  return(list(dataNH1=dataNH1,dataNHim=dataNHim))
}

  dataNH1<-AveIndi_to_PlotLevel(dataNHpi=dataNHpi,dataNHim=dataNHim3yrs_C,popChk=TRUE)$dataNH1
  dataNHim<-AveIndi_to_PlotLevel(dataNHpi=dataNHpi,dataNHim=dataNHim3yrs_C,popChk=TRUE)$dataNHim
#######
#######  
  
######### Only Run for the individual traits
  dataNH<-dataNHim   ### !!!!!  OR dataNH<-dataNH1 ----> In 5.b script RUN
  for (col in morTraits){
    dataNH[,col]<-ifelse(dataNH[,col]==0,NA,dataNH[,col])
  } 
  
### Only run once !!    
  dataNH$bladeLength0<-dataNH$bladeLength
  dataNH$bladeMaxWidth0<-dataNH$bladeMaxWidth
  dataNH$bladeThickness0<-dataNH$bladeThickness
  dataNH$stipeLength0<-dataNH$stipeLength
  dataNH$stipeDiameter0<-dataNH$stipeDiameter
  
  ### if No log
  for (morcol in morTraits){
    dataNH[,morcol]<-log(dataNH[,morcol]+1)
  }
  
  trtlevel<-"IndiLevel"
  traits<-morTraits
  
  dataYr21Plot<-dataYr21
  dataYr21Indi<-droplevels(dataNHim3yrs_C[dataNHim3yrs_C$Year==2021,])
  
  dataYr21<-AveIndi_to_PlotLevel(dataNHpi=dataYr21Plot,dataNHim=dataYr21Indi,popChk=FALSE)$dataNH1
  dataYr21GPlist<-unique(c(as.character(dataYr21$femaPar),as.character(dataYr21$malePar)))
  
### !!!!! Change the corresponding traits for each data set too!!!!
### Run Above 1. vs 2. separatelly, Run only once !!! 
  
  
  
#### Predicting individuals: Add GPs lists (and or SP crosses) as NAs to the phenotypic data
dataNH0<-dataNH
###Crosslist<-rownames(All_grm)[!rownames(All_grm)%in%as.character(dataNH$Crosses)]

# ### If Only Predict the 92 GPs used in Yr2021 SPs !!!!!!!! Results will be the same as predicting a bunch of others
# Crosslist<-dataYr21GPlist
# Trait_Crosses<-unique(as.character(dataNH$Crosses))
# Trait_grm<-All_grm[rownames(All_grm)%in%Trait_Crosses,colnames(All_grm)%in%Trait_Crosses] ### Subset genomic relationship matrix to only include "Crosses" (Yr19+20_SP + Yr21_SP's_GPs) of phenotypic data

### ### If Predict all A_866 individuals  !!!!!!!!
Crosslist<-setdiff(rownames(All_grm),levels(dataNH0$Crosses))
Trait_Crosses<-unique(as.character(dataNH$Crosses))
Trait_grm<-All_grm

dataNAs<-matrix(nrow=length(Crosslist),ncol=ncol(dataNH0))
dataNAs<-as.data.frame(dataNAs)
colnames(dataNAs)<-colnames(dataNH0)
dataNAs$Crosses<-c(Crosslist)
dataNAs$Year<-2019  ### GPs Year NAs used 2019
dataNAs$Year<-as.factor(dataNAs$Year)

dataNH<-rbind(dataNH0,dataNAs)  ### Input phenotypic data Ready.


    ### Input genomic relationship matrix Ready.
    #"SL18-NL-9-FG2xSA18-CB-7-MG1"  "SL18-NC-12-FG2xSL18-NL-5-MG2" in All_grm, Not in Trait_Crosses
  print(dim(dataNH))
  print(sum(is.na(dataNH$Trait)))
  print(sum(rownames(All_grm)%in%Trait_Crosses))
  print(sum(Trait_Crosses%in%rownames(All_grm)))
  print(dim(Trait_grm))   ### ended up being 864,because 2 of them are checks also went into the A_866
  
  ###!!!Later Need to: Add the checks into the Trait_grm, all 1s in diagonal,all 0s for others
  ###"SL18-ME-1-FG1xSL18-ME-1-MG1"  "SL18-ME-1-FG1xSL18-ME-2-MG2" 
  ###"SL18-NL-7-FG3xSL18-OI-15-MG3" "SL18-OI-15-FGxSL18-ME-MG1"
  
  Gencor_Yr<-function(data=data,Trait_grm=Trait_grm4,WithinYr=FALSE){
  
      mod<-asreml(trt~Year+line:Year+block:Year+popChk,
                  random= ~ us(Year):vm(CrossNum,Trait_grm4),
                  residual= ~ dsum(~ idv(units)|Year),
                  data = data, trace=TRUE,
                  maxiter=200,
                  na.action = na.method(y = "include",x="include"))
    ###
    # predMTM = predict(mod, classify = "CrossNum", trace=F)
    ###
    return(list(mod=mod))
  }
  
  modSum<-NULL
  mods<-NULL
  cor<-NULL
  BVs<-NULL # the num row is the length of the bu
  GenOrder<-NULL
  

  for (j in 1:length(traits)){
    
    dataNH$trt<-dataNH[,traits[j]]

    # if (traits[j]%in%c("AshFDwPM")){
    #   dataNH$trt2<-ifelse(dataNH$Year==2019,dataNH$trt*sqrt(10),dataNH$trt)  # Yr1 phenotype * sqrt(10)
    #   dataNH$trt<-dataNH$trt2
    # }

    data<-dataNH    
    ### Adding checks individuals, inside the above dataNH, to the grm
    ChkCross<-unique(droplevels(dataNH[!dataNH$popChk=="ES"&!is.na(dataNH$popChk),])$Crosses) # 33 plots of checks, 4 unique ones
    
    ### Adding the checks diagonal matrix to the G matrix
    Diag<-diag(1,nrow=length(ChkCross),ncol=length(ChkCross))
    Trait_grm4<-magic::adiag(Trait_grm,Diag)
    rownames(Trait_grm4)<-colnames(Trait_grm4)<-c(rownames(Trait_grm),as.character(ChkCross))
    
    
    ### Sort data by year, required by ASReml-R
    data<-data[order(data$Year),]
    
    dataCrosses<-as.character(unique(data$Crosses))
    Trait_grm4<-Trait_grm4[dataCrosses,dataCrosses]  # Order of rows/cols in Grm MUST Match to Phenotypic Crosses level

    ### Create CrossNum, Re-store the Trait_grm4 names to be numeric. Linked in order
    order<-c(1:nrow(Trait_grm4))   
    names(order)<-rownames(Trait_grm4)
    rownames(Trait_grm4)<-colnames(Trait_grm4)<-order
    
    data$CrossNum<-order[as.character(data$Crosses)]
    print(nrow(Trait_grm4)==length(unique(droplevels(data$Crosses))))
    
    ####### Done reordering Trait_grm4
    data$CrossNum<-as.factor(data$CrossNum)
    print(identical(rownames(Trait_grm4),as.character(levels(data$CrossNum)))) #HAS to be TRUE!!!
    
    ####
    mods[[j]]<-Gencor_Yr(data=data,Trait_grm=Trait_grm4,WithinYr=FALSE)$mod
    
    modSum[[j]]<-summary(mods[[j]])$varcomp
    BVs[[j]]<-as.data.frame(mods[[j]]$coefficients$random)  ## Extracting BLUPs/BVs
    BVs[[j]]$order<-rep(order,2)  ## The CrossNum
    BVs[[j]]$CrossOrder<-rep(names(order),2)  ## The actual Crosses
    BVs[[j]]$Year<-substr(rownames(BVs[[j]]),6,9)  ## The Year, all GPs will be Yr2019
    
    cor<-c(cor,modSum[[j]][,"component"][2]/sqrt(modSum[[j]][,"component"][1]*modSum[[j]][,"component"][3]))  # Store Genetic cor for the same trait between Yr19 & Yr20
    GenOrder[[j]]<-order  
  }  ## traits

##### Manually estimate genetic cor, and h2, identical results from vpredict()    
  names(cor)<-traits
    cor
  h2<-matrix(nrow=2,ncol=length(traits))
  colnames(h2)<-traits
  for (j in 1:length(modSum)){
    h2[1,j]<-modSum[[j]][,"component"][1]*2/(modSum[[j]][,"component"][1]*2+modSum[[j]][,"component"][5]) # Yr2019 h2  ##!! The 1,3,5,7 comes from summary(mods[[j]])
    h2[2,j]<-modSum[[j]][,"component"][3]*2/(modSum[[j]][,"component"][3]*2+modSum[[j]][,"component"][7]) # Yr2020 h2  ##!!
  }  
    h2  
  SEh2_Yr19<-NULL  # h2 for Yr19 and se
  SEh2_Yr20<-NULL  # h2 for Yr20 and se
  SEcor<-NULL   #Genetic for for Yr1920 and se
  for (j in 1:length(traits)){
    SEh2_Yr19<-rbind(SEh2_Yr19,vpredict(mods[[j]],h19~V1*2/(V1*2+V5)))    # h2 for Yr2019, vpredict() can automatically generate h2 and its SE from asreml-R
    SEh2_Yr20<-rbind(SEh2_Yr20,vpredict(mods[[j]],h20~V3*2/(V3*2+V7)))    # h2 for Yr2020
    SEcor<-rbind(SEcor,vpredict(mods[[j]],corA~V2/sqrt(V1*V3)))           # Genetic cor for the same trait between years
  }
  
  SingleBV<-matrix(nrow=nrow(BVs[[1]]),ncol=length(traits))  ## convert list into a single data frame
  colnames(SingleBV)<-traits
  rownames(SingleBV)<-stringr::str_split_fixed(rownames(BVs[[1]]),"_",3)[,3]
  for (j in 1:length(traits)){
    SingleBV[,j]<-BVs[[j]]$bu
  }
  corBLUPs<-cor(SingleBV)    # correlation between BLUPs
  SingleBV<-as.data.frame(SingleBV)
  SingleBV$Indiv<-as.vector(BVs[[1]]$CrossOrder)  # The orders in A_866 OR rep(c(names(GenOrder[[1]])),2)
  
  write.csv(cor,paste0(WD,"Genetic_cor_between_Years",trtlevel,".csv"))
  write.csv(SingleBV,paste0(WD,"BVs_UsingYr1920_predict_A866_Indiv_",trtlevel,".csv"))
###### Saved Genetic Cor, h2, se, and the BLUPs
  save(GenOrder,BVs,SEh2_Yr19,SEh2_Yr20,SEcor,SingleBV,mods,file=paste0(WD,"Genetic_cor_h2_between_Years",trtlevel,".Rdata"))

  
###### Load in the BLUEs Yr 2021
  WD<-"/local/workdir/mh865/GCA_SCA/"
  datafdr2<-paste0(WD,"OneTime1920/")
  Yr<-2021
  load(paste0(datafdr2,"BLUEs_Within_",Yr,".Rdata"))
  Trt_BLUE$Crosses<-Trt_BLUE$Row.names
  
  dataYr21<-Trt_BLUE  ### No more chks included
  dataYr21<-droplevels(dataYr21[!duplicated(dataYr21$Crosses),]) ## RM the duplicated crosses, which all have the same BLUPs; 87 ES plotsx27

  corYr21<-matrix(nrow=3,ncol=length(traits))
  AllGEBVs<-matrix(nrow=nrow(BVs[[1]]),ncol=length(traits)) ### GEBVs for Yr19 and Yr20 all traits
  
  Yr21BVfrom19<-matrix(nrow=nrow(dataYr21),ncol=length(traits))
  Yr21BVfrom20<-matrix(nrow=nrow(dataYr21),ncol=length(traits))
  Yr21BVfromComb<-matrix(nrow=nrow(dataYr21),ncol=length(traits))
  rownames(Yr21BVfrom19)<-rownames(Yr21BVfrom20)<-rownames(Yr21BVfromComb)<-dataYr21$Crosses
  
  GPlistPred<-unique(c(as.character(dataYr21$femaPar),as.character(dataYr21$malePar))) ### GPs list for Yr2021
  
  for (j in 1:length(traits)){
    ########### Getting the BVs on GPs
    BV1920<-BVs[[j]]
    BV1920sub<-BV1920[BV1920$CrossOrder%in%GPlistPred,]    ##Subsetting the BLUPs for GPs (Yr21 SP's parents) only
    print(nrow(BV1920sub))
    
    BV1920_wide<-tidyr::spread(BV1920sub,Year,bu)    ## Long to wide formating
    colnames(BV1920_wide)<-c("CrossNum","Crosses","Yr2019","Yr2020") # rename the cols
    
    ### Merging the BVs to the dataYr21 data file
    dataYr21$FGBV19<-expss::vlookup(dataYr21$femaPar,BV1920_wide,lookup_column = "Crosses",result_column = "Yr2019")
    dataYr21$MGBV19<-expss::vlookup(dataYr21$malePar,BV1920_wide,lookup_column = "Crosses",result_column = "Yr2019")
    dataYr21$BV_fromYr19<-dataYr21$FGBV19+dataYr21$MGBV19
    cor19<-cor(dataYr21[,traits[j]],dataYr21$BV_fromYr19,use="complete")
    
    dataYr21$FGBV20<-expss::vlookup(dataYr21$femaPar,BV1920_wide,lookup_column = "Crosses",result_column = "Yr2020")
    dataYr21$MGBV20<-expss::vlookup(dataYr21$malePar,BV1920_wide,lookup_column = "Crosses",result_column = "Yr2020")
    dataYr21$BV_fromYr20<-dataYr21$FGBV20+dataYr21$MGBV20
    
    cor20<-cor(dataYr21[,traits[j]],dataYr21$BV_fromYr20,use="complete")
    
    dataYr21$BVComb1920<-rowMeans(dataYr21[,c("BV_fromYr19","BV_fromYr20")]) ### This is the Index, averaging BVs in Yr2019, Yr2020
    corComb<-cor(dataYr21[,traits[j]],dataYr21$BVComb1920,use="complete")
      
    corYr21[1,j]<-round(cor19,3)
    corYr21[2,j]<-round(cor20,3)
    
    corYr21[3,j]<-round(corComb,3)
    AllGEBVs[,j]<-BVs[[j]]$bu
    
    Yr21BVfrom19[,j]<-dataYr21$BV_fromYr19
    Yr21BVfrom20[,j]<-dataYr21$BV_fromYr20
    Yr21BVfromComb[,j]<-dataYr21$BVComb1920
    
  }
  
  AllGEBVs<-as.data.frame(AllGEBVs)
  
  colnames(AllGEBVs)<-colnames(Yr21BVfrom19)<-colnames(Yr21BVfrom20)<-colnames(Yr21BVfromComb)<-traits
  AllGEBVs$Year<-BVs[[1]]$Year
  AllGEBVs$Crosses<-BVs[[1]]$CrossOrder
  
  colnames(corYr21)<-traits
  rownames(corYr21)<-c("Yr19_Yr21","Yr20_Yr21","MeanTwoYrsBV_Yr21")
  
  write.csv(corYr21,paste0(WD,"GS_cor_Yr1920_predicting_All_A866_Indiv",trtlevel,".csv"))
  save(dataYr21,Yr21BVfrom19,Yr21BVfrom20,Yr21BVfromComb,AllGEBVs,cor,h2,SEh2_Yr19,SEh2_Yr20,SEcor,SingleBV,corYr21,file=paste0(WD,"BVs_UsingYr1920_h2_SE_GenetCor_GScorPredA866Indiv_",trtlevel,".Rdata"))
  
  
#### Cor of BLUPs in different categories of dataYr21. Do it locally on imac
#### load in the BLUPs and BLUEs of Yr2021
  trtlevels<-c("PlotLevel","IndiLevel")
    BV1<-NULL
    BV2<-NULL
    BV3<-NULL
  for (i in 1:length(trtlevels)){
    load(paste0(datafdr,"BVs_UsingYr1920_h2_SE_GenetCor_GScorPredA866Indiv_",trtlevels[i],".Rdata"))  ### This also has the dataYr21
    BV1[[i]]<-Yr21BVfrom19
    BV2[[i]]<-Yr21BVfrom20
    BV3[[i]]<-Yr21BVfromComb
    rm(Yr21BVfrom19)
    rm(Yr21BVfrom20)
    rm(Yr21BVfromComb)
           }
    Yr21BVall1<-merge(BV1[[1]],BV1[[2]],by.x="row.names",by.y="row.names",all.x=TRUE)
    rownames(Yr21BVall1)<-Yr21BVall1$Row.names
    Yr21BVall1<-Yr21BVall1[,!colnames(Yr21BVall1)=="Row.names"]
    
    Yr21BVall2<-merge(BV2[[1]],BV2[[2]],by.x="row.names",by.y="row.names",all.x=TRUE)
    rownames(Yr21BVall2)<-Yr21BVall2$Row.names
    Yr21BVall2<-Yr21BVall2[,!colnames(Yr21BVall2)=="Row.names"]
    
    Yr21BVall3<-merge(BV3[[1]],BV3[[2]],by.x="row.names",by.y="row.names",all.x=TRUE)
    rownames(Yr21BVall3)<-Yr21BVall3$Row.names
    Yr21BVall3<-Yr21BVall3[,!colnames(Yr21BVall3)=="Row.names"]

PlotTraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
morTraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
traits<-c(PlotTraits,morTraits)

CorBLUPs<-function(dataBLUE=dataYr21,dataBV19=Yr21BVall1,dataBV20=Yr21BVall2,dataBVComb=Yr21BVall3){
  cor=matrix(nrow=3,ncol=length(traits))
  for (j in 1:length(traits)){
    cor[1,j]<-cor(dataBLUE[,traits[j]],dataBV19[,traits[j]],use="complete")
    cor[2,j]<-cor(dataBLUE[,traits[j]],dataBV20[,traits[j]],use="complete")
    cor[3,j]<-cor(dataBLUE[,traits[j]],dataBVComb[,traits[j]],use="complete")
  }
  colnames(cor)<-traits
  rownames(cor)<-c("BVfromYr19","BVfromYr20","BVfromComb")
  return(cor)
}  

### Make order the same!!
dataYr21<-dataYr21[match(rownames(Yr21BVall1),dataYr21$Row.names),] 
identical(as.character(rownames(Yr21BVall1)),as.character(dataYr21$Row.names))   ### !!! Has to be TRUE
CorBLUPs(dataBLUE=dataYr21,dataBV19=Yr21BVall1,dataBV20=Yr21BVall2,dataBVComb=Yr21BVall3)


### Load in the original data where there is "category"
WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003"
datafdr<-paste0(WD,"/data/")
Yr2021<-read.csv(paste0(datafdr,"2021_PhenotypingDatasheet_PlotLevel.csv"),sep=",",header=TRUE)
#RM other crosses categories
RMCat<-c("","1000","2000","250","4000","500","Industry Cross","Paint_replace","Spray","Spray_replace")
Yr2021_sub1<-droplevels(Yr2021[!Yr2021$Category%in%RMCat,colnames(Yr2021)%in%c("Crosses","Category","Photo.Score","Cross.ID")])   
  dim(Yr2021_sub1)  #280

dataYr21$Category<-expss::vlookup(dataYr21$Crosses,Yr2021_sub1,result_column = "Category",lookup_column = "Crosses")

table(droplevels(dataYr21$Category))

#SubX<-dataYr21[dataYr21$Category=="Top160",]$Crosses   ### Predict the top 160X

SubX<-dataYr21[dataYr21$Category=="Extra",]$Crosses
identical(as.character(dataYr21$Row.names),as.character(rownames(Yr21BVall1)))  ### !!! Has to be TRUE
CorBLUPs(dataBLUE=dataYr21[dataYr21$Row.names%in%SubX,],dataBV19=Yr21BVall1[rownames(Yr21BVall1)%in%SubX,],dataBV20=Yr21BVall2[rownames(Yr21BVall2)%in%SubX,],dataBVComb=Yr21BVall3[rownames(Yr21BVall3)%in%SubX,])   
 
### !!!!!!!!!
### Load in the BLUEs of Yr2020, find those plots that are common betwen Yr21 and Yr20
Yr<-2020
WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003"
datafdr<-paste0(WD,"/data/")
load(paste0(datafdr,"BLUEs_Within_",Yr,".Rdata"))
Trt_BLUE$Crosses<-Trt_BLUE$Row.names
dataYr20<-Trt_BLUE

ComX<-dataYr20$Row.names[dataYr20$Row.names%in%dataYr21$Row.names]  ### Correlating the 12 common plots

ComX20<-dataYr20[dataYr20$Row.names%in%ComX,]
ComX21<-dataYr21[dataYr21$Row.names%in%ComX,]
identical(as.character(ComX20$Row.names),as.character(ComX21$Row.names))
identical(as.character(ComX21$Row.names),as.character(rownames(Yr21BVall1)[rownames(Yr21BVall1)%in%ComX]))
CorBLUPs(dataBLUE=ComX21,dataBV19=Yr21BVall1[rownames(Yr21BVall1)%in%ComX,],dataBV20=Yr21BVall2[rownames(Yr21BVall2)%in%ComX,],dataBVComb=Yr21BVall3[rownames(Yr21BVall3)%in%ComX,])   


cor<-matrix(nrow=1,ncol=length(traits))
for (j in 1:length(traits)){
  cor[1,j]<-cor(ComX20[,traits[j]],ComX21[,traits[j]],use="complete")
}
  colnames(cor)<-traits
  cor


  ########### Compare the BVs, when GPs NAs Year used Yr2019 OR Yr2020 
  # BVs<-BVs[[j]]
  # BVs$order<-rep(order,2)
  # BVs$CrossOrder<-rep(names(order),2)
  # write.csv(BVs,"CrossList870_UseYr2019forNAs_pred_YR21.csv")
  #  BV19<-read.csv("CrossList870_UseYr2019forNAs_pred_YR21.csv",sep=",",header=TRUE)
  #  BV20<-read.csv("CrossList870_UseYr2020forNAs_pred_YR21.csv",sep=",",header=TRUE)
  # 
  #  BV19sub<-BV19[BV19$CrossOrder%in%Crosslist,]
  #  BV20sub<-BV20[BV20$CrossOrder%in%Crosslist,]
  #  
  #  BV19sub$Year<-substr(BV19sub$X,6,9)
  #  BV20sub$Year<-substr(BV20sub$X,6,9)
  #  
  #  BV19sub$ID<-paste0(BV19sub$Year,"_",BV19sub$CrossOrder)
  #  BV20sub$ID<-paste0(BV20sub$Year,"_",BV20sub$CrossOrder)
  #  
  #  BV1920<-merge(BV19sub,BV20sub,by.x="ID",by.y="ID",all.x=TRUE)
  #  cor(BV1920$b8.x,round(BV1920$bu.y,2))  #1
  # 
  #  #### Long to Wide
  # # library(tidyr)
  #  BV1920sub<-BV1920[,c("bu.x","order.x","CrossOrder.x","Year.x")]
  #  BV1920_wide<-spread(BV1920sub,Year.x,bu.x) 
  #  colnames(BV1920_wide)<-c("Order_Yrused2019","CrossOrder","Yr2019","Yr2020")
  #  cor(BV1920_wide$Yr2019,BV1920_wide$Yr2020)  #0.977
  
  ###########
  # This below piece Did not work to average across years
  # BVCombined<-predict(mod,classify="CrossNum",trace=T) ## NA used Yr2019
  # BVCombined2<-predict(mod2,classify="CrossNum",trace=T) ## NA used Yr2020