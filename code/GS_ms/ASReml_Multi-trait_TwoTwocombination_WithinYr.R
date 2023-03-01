#### ASReml multi-trait genetic cor, Two two combination wihin each year analysis. Using BLUE !!! ----- Did not use in the ms in the end
#### Yr 2021, perDW and Ash 2-2 comb, Var structure not positive definite

### Results Output are in:"/local/workdir/mh865/GCA_SCA/OneTime1920/data/"

rm(list=ls())

#### Data sort by Year
#### Read in the A matrix data
#### reverse the A matrix
#### Use model

WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal
datafdr<-paste0(WD,"OneTime1920/data/")

PlotTraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades")  
morTraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
traits<-c(PlotTraits,morTraits)

#### 1.Read in the BLUE data
datafdr2<-paste0(WD,"OneTime1920/")
AllBLUEs<-NULL
for (Yr in c(2019,2020,2021)){
  
  load(file=paste0(datafdr2,"BLUEs_Within_",Yr,".Rdata"))
  Y<-droplevels(Trt_BLUE)  ### droplevels so that the # of crosses were matching!
  Y$Crosses<-as.factor(Y$Row.names)
  Y<-Y[!duplicated(Y$Crosses),]
  
  AllBLUEs<-rbind(AllBLUEs,Y)
  }
AllBLUEs$Year<-as.factor(AllBLUEs$Year.x)
#### BLUEs were already log transformed

#### 2.Load the A matrix, outCovComb
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021_866Indiv.Rdata")) ##!!!
A_866<-outCovComb4_MixOrder
#load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0712_2021.Rdata")) ##!!! 950x950 Indiv


Plottraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
Mortraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
Alltraits<-c(Plottraits,Mortraits)

TraitComb<-NULL   ### Traits two-two combinations
for (i in 1:(length(Alltraits)-1)){
  trt1<-Alltraits[i]
  trt2<-Alltraits[c((i+1):length(Alltraits))]
  trtMatch<-rep(trt1,length(trt2))
  dftmp<-cbind(trtMatch,trt2)
  TraitComb<-rbind(TraitComb,dftmp)
}
TraitComb<-as.data.frame(TraitComb)
TraitComb$cor<-rep(NA,nrow(TraitComb))


for (Yr in c(2019, 2020)){
  dataNHpi<-droplevels(AllBLUEs[AllBLUEs$Year==Yr,])
  dataNH<-dataNHpi
  Trait_Crosses<-as.character(dataNH$Crosses)
  Trait_grm<-A_866[rownames(A_866)%in%Trait_Crosses,colnames(A_866)%in%Trait_Crosses]
    print(dim(dataNH))
    print(sum(is.na(dataNH$Trait)))
    print(dim(Trait_grm))
   
### Convert Crosses to CrossNum
  data<-dataNH
    print(nrow(Trait_grm)==length(unique(droplevels(data$Crosses))))
    print(identical(rownames(Trait_grm),as.character(levels(data$Crosses))))
    
    Cross_link<-as.data.frame(rownames(Trait_grm))
    colnames(Cross_link)<-"Cross_chr"
    Cross_link$num<-1:nrow(Trait_grm)
    
    ### Look up the CrossNum for A
    rownames(Trait_grm)<-expss::vlookup(rownames(Trait_grm),dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )
    colnames(Trait_grm)<-expss::vlookup(colnames(Trait_grm),dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )
    
    ### look up the CrossNum for data
    data$CrossesNum<-expss::vlookup(data$Crosses,dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )
    data$CrossesNum<-as.factor(data$CrossesNum)
    
    ### Order the data$CrossesNum to be matching with that in Trait_grm3 CrossesNum 
    data<-data[order(data$Year,data$CrossesNum),]
    
    print(identical(rownames(Trait_grm),as.character(levels(data$CrossNum))))
    
    ### Setting up models
    mods<-NULL 
  for (k in 1:nrow(TraitComb)) {
    data$trt1<-data[,TraitComb[k,1]]
    data$trt2<-data[,TraitComb[k,2]]

    ### Using Only 2 traits
    modMT1<-asreml(cbind(trt1,trt2) ~ trait,
                   random= ~ us(trait):vm(CrossesNum,Trait_grm),
                   residual = ~id(units):us(trait),
                   data = data, maxiter=100, trace=TRUE)
    
    GenVar1<-summary(modMT1)$varcomp ###!!! Codes below can be used for different models
    #write.csv(GenVar1,"Plot5traits_Indivial_MultiTrait_Genetic_var_04222021.csv") #and_
    
    colnames(GenVar1)
    ## Extract rows with extra parts
    
    Gen_varcov<-GenVar1[grepl("trait:vm",rownames(GenVar1)),]
    ## RM names with extra parts
    rownames(Gen_varcov)<-stringr::str_split_fixed(rownames(Gen_varcov),"_",3)[,3]
    
    #Mtraits<-c("wetWgtPlot","dryWgtPerM","AshFDwPM","densityBlades")
    
    ### The varcomp is in long format, convert to wide format
    Gen_varcov$var1<-stringr::str_split_fixed(rownames(Gen_varcov),":",2)[,1]
    Gen_varcov$var2<-stringr::str_split_fixed(rownames(Gen_varcov),":",2)[,2]
    Gen_varcov$cor<-Gen_varcov$component
    Gen_varcov<-Gen_varcov[order(rownames(Gen_varcov)),c("var1","var2","component")]
    Gen_varcov_wide<-tidyr::spread(Gen_varcov,var2,component)
    
    ### Reorganize: Fill in NAs
    rownames(Gen_varcov_wide)<-Gen_varcov_wide$var1
    Gen_varcov_wide<-Gen_varcov_wide[,-1]
    
    Gen_varcov_imp<-matrix(nrow=nrow(Gen_varcov_wide),ncol=ncol(Gen_varcov_wide))
    rownames(Gen_varcov_imp)<-rownames(Gen_varcov_wide)
    colnames(Gen_varcov_imp)<-colnames(Gen_varcov_wide)
    
    for (i in 1:nrow(Gen_varcov_wide)){
      for (j in 1:ncol(Gen_varcov_wide)){
        Gen_varcov_imp[i,j]<-ifelse(is.na(Gen_varcov_wide[i,j]),Gen_varcov_wide[j,i],Gen_varcov_wide[i,j])
      }
    }
    
    print(isSymmetric(round(Gen_varcov_imp,5)))
    Gen.cor<-cov2cor(Gen_varcov_imp)
    
    TraitComb$cor[k]<-Gen.cor[1,2]
    mods[[k]]<-modMT1
  } ### TraitComb two-two combinations
    
    
    #### Using more than 2 traits
    modMT2<-asreml(cbind(dryWgtPerM,bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait,
                   random= ~ us(trait):vm(CrossesNum,Trait_grm),
                   residual = ~id(units):us(trait),
                   data = data, maxiter=100, trace=TRUE)
    
    
    
    save(mods,file=paste0(datafdr,"Alltraits_Genetic_cor_Multi-variate_BLUE_Yr",Yr,"_TwoTwoComb_modOutput.Rdata"))
    Genetic_cor<-tidyr::spread(TraitComb,trt2,cor)
    
    rownames(Genetic_cor)<-Genetic_cor$trtMatch
    Genetic_cor<-Genetic_cor[,!colnames(Genetic_cor)=="trtMatch"]
    
    ordersRow<-Alltraits[-length(Alltraits)]  # No last variable in the rownames of Genetic_cor
    ordersCol<-Alltraits[-1]  # No first variable in the colnames of Genetic_cor
    Genetic.cor<-Genetic_cor[match(ordersRow,rownames(Genetic_cor)),match(ordersCol,colnames(Genetic_cor))]

    write.csv(Genetic.cor,paste0(datafdr,"Alltraits_Genetic_cor_Multi-variate_BLUE_Yr",Yr,".csv")) #
}



#### Compare results from BLUE vs Raw data modeling
Yr<-2019
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal
datafdr<-paste0(WD,"OneTime1920/data/")

BLUEMulticor<-read.csv(paste0(datafdr,"Alltraits_Genetic_cor_Multi-variate_BLUE_Yr",Yr,".csv"))
RawMulticor<-read.csv(paste0(datafdr,"Alltraits_Genetic_cor_Multi-variate_RawData_Yr",Yr,".csv"))

identical(colnames(BLUEMulticor),colnames(RawMulticor))
identical(rownames(BLUEMulticor),rownames(RawMulticor))
cor(c(as.matrix(BLUEMulticor[,-1])),c(as.matrix(RawMulticor[,-1])),use="complete")   ### cor ~ 0.91

# ### Inverse the Trait_grm
# snpRelMat<-Trait_grm
# Gsnp=solve(snpRelMat+diag(1e-6,length(snpRelMat[,1]),length(snpRelMat[,1])))
# 
# nrow(Gsnp)==sum(rownames(Gsnp)==levels(data$CrossesNum))
# nrow(Gsnp)==sum(colnames(Gsnp)==levels(data$CrossesNum))
# 
# attr(Gsnp, "INVERSE")=TRUE
# asreml.options(extra=30)