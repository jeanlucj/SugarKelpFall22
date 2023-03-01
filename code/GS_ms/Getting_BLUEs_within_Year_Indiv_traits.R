### II. Getting BLUEs within Year, now for Indiv traits
### III. Merged the BLUEs from Plot level to Indi traits

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-paste0(WD,"OneTime192021/data/")
load(paste0(datafdr,"3yrs_data_plotformated_indivNot_12012021.Rdata"))    ### Included the Yr2021 SPs in pedigree
datafdr2<-paste0(WD,"OneTime1920/")
dataNHim<-dataNHim3yrs_C

morTraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

dataNHim$line<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "line")
dataNHim$block<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "block")
dataNHim$Year<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Year")
dataNHim$popChk<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "popChk")
dataNHim$Crosses<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Crosses")
dataNHim$PhotoScore<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "PhotoScore")
dataNHim<-dataNHim[which(dataNHim$PhotoScore >1),]  # 3969 rows with PhotoScore >1

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

### log transform
for (morcol in morTraits){
  dataNH[,morcol]<-log(dataNH[,morcol]+1)
}

trtlevel<-"IndiLevel"
traits<-morTraits
### Only run once !!


### Creating the group factor, that separates the Checks and ES plots
dataNH$group<-ifelse(dataNH$popChk=="ES",1,0) # 1 is the ES
dataNH$group<-as.factor(dataNH$group)

data<-dataNH

library(stringr)
library(asreml)
library(ggplot2)  
library(expss)

for (i in 1:length(traits)){
  data$Trait<-data[,traits[i]]
  dev.set()
  pdf(paste0(traits[i],"_threeyears.pdf"))
  plot2<-ggplot(data=data,aes(Trait,Year))+
    geom_point(aes(color=as.factor(Year)))+ 
    geom_line(aes(group=as.factor(Crosses)))
  print(plot2)
  dev.off()
}

for (Yr in c(2019, 2020, 2021)){
  dataYr0<-droplevels(data[data$Year==Yr,])
 
  dataYr<-dataYr0
  corBLUEBLUPs<-NULL
  BLUEdf<-matrix(nrow=nlevels(dataYr$Crosses),ncol=length(traits))
  rownames(BLUEdf)<-as.character(unique(dataYr$Crosses))
  colnames(BLUEdf)<-traits
  
  deBLUPdf<-matrix(nrow=nlevels(dataYr$Crosses),ncol=length(traits))
  rownames(deBLUPdf)<-as.character(unique(dataYr$Crosses))
  colnames(deBLUPdf)<-traits
  mod_P<-NULL
  mod_E<-NULL
  corEP<-NULL
for (i in 1:length(traits)){ 
  
  dataYr$Trait<-dataYr[,traits[i]]
  
  if (Yr %in% c(2019,2020)){
    mod3<-asreml(Trait~ popChk+line+block+Crosses,
                 residual= ~ idv(units),
                 data = dataYr, trace=TRUE,
                 maxiter=200,
                 na.action = na.method(y = "include",x="include"))
    
    mod4<-asreml(Trait~ popChk+line+block,
                 random= ~ Crosses,
                 residual= ~ idv(units),
                 data = dataYr, trace=TRUE,
                 maxiter=200,
                 na.action = na.method(y = "include",x="include"))
  }else if (Yr ==2021){
    mod3<-asreml(Trait~ line+block+Crosses,
                 residual= ~ idv(units),
                 data = dataYr, trace=TRUE,
                 maxiter=200,
                 na.action = na.method(y = "include",x="include"))
    
    mod4<-asreml(Trait~ line+block,
                 random= ~Crosses,
                 residual= ~ idv(units),
                 data = dataYr, trace=TRUE,
                 maxiter=200,
                 na.action = na.method(y = "include",x="include"))
    
  }
  
  ### Getting BLUE
  mod_E[[i]]<-mod3
  
  modE<-mod3
  BLUE<-as.data.frame(summary(modE,coef=TRUE)$coef.fixed)
  BLUE<-BLUE[grepl("Crosses_",rownames(BLUE)),]
  BLUE$Crosses<-str_split_fixed(rownames(BLUE),"_",2)[,2]
  ## cor(mod$coefficients$fixed,as.data.frame(summary(mod,coef=TRUE)$coef.fixed)$solution) #1
  
  BLUEdf[,i]<-expss::vlookup(rownames(BLUEdf),dict=BLUE,lookup_column="Crosses",result_column = "solution")
  
  print(Yr)
  dev.set()
  pdf(paste0("EvaluateMod3_BLUE","_",traits[i],"_",Yr,".pdf"))
  plot(mod3)
  dev.off()
  
  dev.set()
  pdf(paste0("EvaluateMod4_BLUP","_",traits[i],"_",Yr,".pdf"))
  plot(mod4)
  dev.off()
  
  dev.set()
  pdf(paste0(traits[i],"_",Yr,"_hist.pdf"))
  plot(hist(dataYr$Trait))
  dev.off()
  
  ### Getting deBLUPs
  mod_P[[i]]<-mod4
  
  modP<-mod4
  BLUPs<-as.data.frame(summary(modP,coef=TRUE)$coef.random)
  BLUPs<-BLUPs[grepl("Crosses_",rownames(BLUPs)),]
  BLUPs$Crosses<-str_split_fixed(rownames(BLUPs),"_",2)[,2]
  BLUPs$PEV<-BLUPs[,"std.error"]^2
  BLUPs$TotalGV<-summary(modP)$varcomp["Crosses","component"]
  BLUPs$tmpBLUPs<-BLUPs$solution  
  BLUPs$dBLUPs<-BLUPs$tmpBLUPs/(1-BLUPs$PEV/BLUPs$TotalGV)
  
  ### Comparing BLUE VS deBLUPs
  deBLUPdf[,i]<-expss::vlookup(rownames(deBLUPdf),dict=BLUPs,lookup_column="Crosses",result_column = "dBLUPs")
    print(identical(BLUE$Crosses,BLUPs$Crosses))
  corBLUEBLUPs<-c(corBLUEBLUPs,cor(BLUE$solution,BLUPs$dBLUPs)) 
  
  deBLUPES<-deBLUPdf[rownames(deBLUPdf)%in%dataYr[dataYr$popChk=="ES",]$Crosses,]
  BLUEdfES<-BLUEdf[rownames(BLUEdf)%in%dataYr[dataYr$popChk=="ES",]$Crosses,]
  corEP<-c(corEP,cor(BLUEdfES[,i],deBLUPES[,i],use="complete"))
  
      } 
  print(corEP)
  save(BLUEdf,BLUEdfES,deBLUPdf,deBLUPES,mod_E,mod_P,corEP,dataYr,file=paste0("deBLUP_BLUE_",trtlevel,"_",Yr,".Rdata"))
  
}

#### Use wald(mod_E[[i]]) to check the significance of the models variables


#### III. Combine both Plotlevel and Indilevel traits
colkeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","group")

PlotTraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades")  
morTraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
Alltraits<-c(PlotTraits,morTraits)

for (Yr in c(2019,2020,2021)){
    
### Load Indi level traits, Deal with NAs  
  load(paste0("deBLUP_BLUE_","IndiLevel","_",Yr,".Rdata"))
  for (i in 1:length(morTraits)){
    NACrosses<-dataYr$Crosses[is.na(dataYr[,morTraits[i]])] 
    BLUEdfES[,morTraits[i]]<-ifelse(rownames(deBLUPES)%in%NACrosses,NA,deBLUPES[,morTraits[i]])
  } ### If Original data is NA, then BLUEs also asign to be NAs
  
  Indi_BLUE<-BLUEdfES  ## Only took the ES plots, for Indi level
  rm(BLUEdfES)
  print(Yr)
  print(corEP)
  
### Load Plot and Indi level BLUE, Deal with NAs    
  load(paste0("deBLUP_BLUE_","PlotLevel","_",Yr,".Rdata"))
  for (i in 1:length(PlotTraits)){
    NACrosses<-dataYr$Crosses[is.na(dataYr[,PlotTraits[i]])]
    BLUEdfES[,PlotTraits[i]]<-ifelse(rownames(deBLUPES)%in%NACrosses,NA,deBLUPES[,PlotTraits[i]])
  } ### If Original data is NA, then BLUEs also asign to be NAs
  
  Plot_BLUE<-BLUEdfES  ## Only took the ES plots, for Plot level /// CHANGE if use deBLUPs
  rm(BLUEdfES)
  print(corEP)
  
### Merge the BLUE of Plot and Indi level  
  Exp_fct<-dataYr[,colkeep]  ## Experimental factors from the Plot_level dataYr
  print(dim(dataYr));print(nlevels(dataYr$Crosses))
  
  Both_BLUE<-merge(Plot_BLUE,Indi_BLUE,by.x="row.names",by.y="row.names",all.x=TRUE)
  Both_BLUE$Year<-Yr
  Trt_BLUE<-merge(Both_BLUE,Exp_fct,by.x="Row.names",by.y="Crosses",all.x=TRUE) ##??
  
  save(Trt_BLUE,file=paste0(datafdr2,"BLUEs_Within_",Yr,".Rdata"))  ## BLUEs for both Plot and Indi
}

#### Use wald(mod_E[[i]]) to check the significance of the models variables



# ### compare sommer deBLUP to asreml deBLUP 
# load("/local/workdir/mh865/GCA_SCA/OneTime1920/data/Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD_08182021.Rdata")
# Yr<-2020
# load(paste0("/local/workdir/mh865/ASreml/","deBLUP_BLUE_",trtlevel,"_",Yr,".Rdata"))
# 
# SommerdeBLUPs<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,]
# SommerdeBLUPs<-SommerdeBLUPs[match(rownames(deBLUPES),SommerdeBLUPs$Crosses.x),]
# identical(SommerdeBLUPs$Crosses.x,rownames(deBLUPES))
# 
# cors<-NULL
# for (i in 1:length(traits)){
#   cors<-c(cors,cor(SommerdeBLUPs[,traits[i]],deBLUPES[,i],use="complete"))
# }
# # cors: 0.9525276 0.9470991 0.9909628 0.8976982 0.9553857
# # corEP:
##2019
#corBLUE_deBLUP Plot 0.9203303 0.7997792 0.7389927 0.8951899 0.9458114
#corBLUE_deBLUP Indi 0.9199891 0.7356590 0.7450783 0.9899055 0.8334329 0.8691223

## 2020
#corBLUE_deBLUP Plot 0.8773258 0.8065603 0.8051721 0.8655930 0.7981260
#corBLUE_deBLUP indi 0.9251436 0.8244757 0.9560670 0.5563315 0.8477388 0.7081772

##2021
#corBLUE_deBLUP Plot 0.6001929 0.4818751 0.6880880 0.8689542 0.9061190
#corBLUE_deBLUP Indi 0.6756775 0.4005580 0.5552061 0.3111085 0.7143073 0.8057074



# ### If just use 90% data:
# a<-round(0.6*nrow(dataYr0))
# dataYr<-droplevels(dataYr0[c(1:a,(nrow(dataYr0)-18):nrow(dataYr0)),])

# ### Comparing DIFFERENT models
# mod3.5<-asreml(Trait~ line+Crosses,
#                random=~block,
#                residual= ~ idv(units),
#                data = dataYr, trace=TRUE,
#                maxiter=200,
#                na.action = na.method(y = "include",x="include"))
# mod3.6<-asreml(Trait~ Crosses,
#                random=~line+block,
#                residual= ~ idv(units),
#                data = dataYr, trace=TRUE,
#                maxiter=200,
#                na.action = na.method(y = "include",x="include"))  # mod3.6 similar AIC,BIC to mod3.5; mod3 better BIC/AIC than mod3.5
# 
# mod4.5<-asreml(Trait~ line,
#                random= ~Crosses+block,
#                residual= ~ idv(units),
#                data = dataYr, trace=TRUE,
#                maxiter=200,
#                na.action = na.method(y = "include",x="include"))
# mod4.6<-asreml(Trait~ 1,
#                random= ~Crosses+block+line,
#                residual= ~ idv(units),
#                data = dataYr, trace=TRUE,
#                maxiter=200,
#                na.action = na.method(y = "include",x="include"))  # block not estimable
