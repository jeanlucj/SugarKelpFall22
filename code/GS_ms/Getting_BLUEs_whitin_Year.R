#### I. Getting deBLUP/BLUEs from within Year for Plot Traits
#### load in the 3yrs data for plot level
#### Either deBLUP or BLUEs to GS. deBLUPs should be weighed. I changed to use BLUEs. It would also be easy use BLUEs directly as any phenotypic value comparisons
### dataNHpi3yrs_C_Ash is the raw 3yrs data without popChk
### dataNHpi_RMchk is the raw 3yrs data removed checks

### The analysis Use these two:
### dataNHpi is the formated raw 3yrs data with all experimental factors
### dataNHim3yrs_C is the non-formated raw 3yrs data for Individual level traits ONLY
### Yr2019 and Yr2021: Ash/AshFDWpM has to remove line/blocks in the model. Yr2021: All traits have to remove popChk in the model.

#### load in the 3yrs data for individual level
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-paste0(WD,"OneTime192021/data/")
load(paste0(datafdr,"3yrs_data_plotformated_indivNot_12012021.Rdata"))    ### Included the Yr2021 SPs in pedigree
dataNH<-dataNHpi   # dataNHpi is the 3 yrs plot level data already formated

PlotTraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades")  


### Only Run ONCE!!!!
dataNH$wetWgtPerM0<-dataNH$wetWgtPerM
dataNH$percDryWgt0<-dataNH$percDryWgt
dataNH$dryWgtPerM0<-dataNH$dryWgtPerM
dataNH$Ash0<-dataNH$Ash
dataNH$AshFDwPM0<-dataNH$AshFDwPM
dataNH$densityBlades0<-dataNH$densityBlades

### log transform
for (col in PlotTraits){
  dataNH[,col]<-log(dataNH[,col]+1)
}

trtlevel<-"PlotLevel"
traits<-PlotTraits
### Only Run ONCE!!!!


### Creating the group factor, that separates the Checks and ES plots
dataNH$group<-ifelse(dataNH$popChk=="ES",1,0) # 1 is the ES
dataNH$group<-as.factor(dataNH$group)

dataNH$line<-as.factor(dataNH$line)
dataNH$block<-as.factor(dataNH$block)

data<-dataNH

library(stringr)
library(asreml)
library(ggplot2)  

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
    
    ## If just use 90% data:
    #a<-round(0.6*nrow(dataYr0))
    #dataYr<-droplevels(dataYr0[c(1:a,(nrow(dataYr0)-18):nrow(dataYr0)),])
    
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
    
 if (Yr==2019){
    ### mod3:Getting BLUEs mod4:Getting BLUPs
       if (traits[i]%in% c("Ash","AshFDwPM")){
      
      mod3<-asreml(Trait~ popChk+Crosses,
                     residual= ~ idv(units),
                     data = dataYr, trace=TRUE,
                     maxiter=200,
                     na.action = na.method(y = "include",x="include"))
      
      mod4<-asreml(Trait~ popChk,
                   random= ~Crosses,
                   residual= ~ idv(units),
                   data = dataYr, trace=TRUE,
                   maxiter=200,
                   na.action = na.method(y = "include",x="include"))
      
        }else{
      mod3<-asreml(Trait~ popChk+line+block+Crosses,
                 residual= ~ idv(units),
                 data = dataYr, trace=TRUE,
                 maxiter=200,
                 na.action = na.method(y = "include",x="include"))
      
      mod4<-asreml(Trait~ popChk+line+block,
                   random= ~Crosses,
                   residual= ~ idv(units),
                   data = dataYr, trace=TRUE,
                   maxiter=200,
                   na.action = na.method(y = "include",x="include"))
      ### For percDryWgt in Yr2019, popChk works; but popChk gives a "B" (varcomp summary not bounding) for the Crosses variable 
      
           }
  }else if (Yr==2020){
    
      mod3<-asreml(Trait~ popChk+line+block+Crosses,
                   residual= ~ idv(units),
                   data = dataYr, trace=TRUE,
                   maxiter=200,
                   na.action = na.method(y = "include",x="include"))
         
      mod4<-asreml(Trait~ popChk+line+block,
                   random= ~Crosses,
                   residual= ~ idv(units),
                   data = dataYr, trace=TRUE,
                   maxiter=200,
                   na.action = na.method(y = "include",x="include"))
      
                    
  }else if (Yr==2021){
        # No popchk
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
    
    ### Getting the BLUEs
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
    
    ### Getting the deBLUPs
    mod_P[[i]]<-mod4
    
    modP<-mod4
    
    BLUPs<-as.data.frame(summary(modP,coef=TRUE)$coef.random)
    BLUPs<-BLUPs[grepl("Crosses_",rownames(BLUPs)),]
    BLUPs$Crosses<-str_split_fixed(rownames(BLUPs),"_",2)[,2]
    BLUPs$PEV<-BLUPs[,"std.error"]^2
    BLUPs$TotalGV<-summary(modP)$varcomp["Crosses","component"]
    BLUPs$tmpBLUPs<-BLUPs$solution  
    BLUPs$dBLUPs<-BLUPs$tmpBLUPs/(1-BLUPs$PEV/BLUPs$TotalGV)
    
    ### Compare BLUE and the deBLUPs
    deBLUPdf[,i]<-expss::vlookup(rownames(deBLUPdf),dict=BLUPs,lookup_column="Crosses",result_column = "dBLUPs")
    identical(BLUE$Crosses,BLUPs$Crosses) ### Must be true
    corBLUEBLUPs<-c(corBLUEBLUPs,cor(BLUE$solution,BLUPs$dBLUPs))  # 0.87 :WWpM
    
    ### Subsetting only the ES plots
    deBLUPES<-deBLUPdf[rownames(deBLUPdf)%in%dataYr[dataYr$popChk=="ES",]$Crosses,]
    BLUEdfES<-BLUEdf[rownames(BLUEdf)%in%dataYr[dataYr$popChk=="ES",]$Crosses,]
    corEP<-c(corEP,cor(BLUEdfES[,i],deBLUPES[,i],use="complete"))
} 
    print(corEP)
    ### Raw data NAs will be assigned to 0 in BLUEs. Will convert later
    save(BLUEdf,BLUEdfES,deBLUPdf,deBLUPES,mod_E,mod_P,corEP,dataYr,file=paste0("deBLUP_BLUE_",trtlevel,"_",Yr,".Rdata"))
}

#### Use wald(mod_E[[i]]) to check the significance of the models variables



# #### Compare the deBLUPs to those from sommer for Yr2019
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


###### A different way of getting the BLUPs
## BLUPs2<-predict(mod,classify="Crosses")[[1]]
## identical(BLUPs$Crosses,as.character(BLUPs2$Crosses))
## cor(BLUPs[,"solution"],BLUPs2$predicted.value)  ## 1

###### Run model with only 60% of the data
# library(lme4)
# library(lmerTest)
#  for (i in 1:length(traits)){
#    dataYr$Trait<-dataYr[,traits[i]]
#    modTest<-lmer(Trait~group+(1|Crosses)+(1|line),data=dataYr)
#  }
# ranova(modTest) # Testing random # singularity issues for line and blocks within Yr2021 (no enough data points to estimate)

## BLUE and deBLUPs in ASReml corBLUEBLUPs (included check crosses)
## 60% data corBLUEBLUPs: 0.8991823 0.9505846 0.9420331 0.9494987 0.8361907 0.9069070
## 100 % data corBLUEBLUPs: 0.9230970 0.9518619 0.9498605 0.9898118 0.8322359 0.9187706


#### DIFFERENT Model comparisons ---- Picked model 3.8
#### Compared to mod3
# dataYr$Trait<-dataYr[,traits[i]]
# 
mod13<-asreml(fixed=Trait~at(group,'1'):Crosses+at(group,'0'):Crosses,
                  random= ~line+block,
                 residual=~idv(units),data=dataYr,na.action = na.method(y = "include",x="include"))
summary(mod13)

mod3<-asreml(Trait~ popChk+Crosses,
             residual= ~ idv(units),
             data = dataYr, trace=TRUE,
             maxiter=200,
             na.action = na.method(y = "include",x="include"))

mod3.5<-asreml(Trait~ group+Crosses,
               residual= ~ idv(units),
               data = dataYr, trace=TRUE,
               maxiter=200,
               na.action = na.method(y = "include",x="include"))

mod3.6<-asreml(Trait~ group+Crosses,
               random= ~ line+block,
               residual= ~ idv(units),
               data = dataYr, trace=TRUE,
               maxiter=200,
               na.action = na.method(y = "include",x="include"))

mod3.7<-asreml(Trait~ group+Crosses+line+block,
               residual= ~ idv(units),
               data = dataYr, trace=TRUE,
               maxiter=200,
               na.action = na.method(y = "include",x="include"))

mod3.8<-asreml(Trait~ popChk+Crosses+line+block,
               residual= ~ idv(units),
               data = dataYr, trace=TRUE,
               maxiter=200,
               na.action = na.method(y = "include",x="include"))

BLUEs<-function(mod){
  BLUE<-as.data.frame(summary(mod,coef=TRUE)$coef.fixed)
  BLUE<-BLUE[grepl("Crosses_",rownames(BLUE)),]
  BLUE$Crosses<-str_split_fixed(rownames(BLUE),"_",2)[,2]
  
  NACrosses<-dataYr$Crosses[is.na(dataYr[,traits[i]])] 
  BLUE$BLUE<-ifelse(BLUE$Crosses%in%NACrosses,NA,BLUE$solution)
  
  return(BLUE)
}

BLUE3<-BLUEs(mod=mod3)

BLUE3.5<-BLUEs(mod=mod3.5)

BLUE3.6<-BLUEs(mod=mod3.6)

BLUE3.7<-BLUEs(mod=mod3.7)

BLUE3.8<-BLUEs(mod=mod3.8)

mod3.9<-lm(Trait~ popChk+Crosses+line+block,data=dataYr) # compare to mod3.8
BLUE3.9<-as.data.frame(lsmeans(mod3.9,"Crosses"))
BLUE3.9$BLUE<-BLUE3.9$lsmean

mod3.10<-lmer(Trait~popChk+Crosses+(1|line),data=dataYr)
BLUE3.10<-as.data.frame(lsmeans(mod3.10,"Crosses"))

cors<-function(data1,data2){
  data1$BLUE2<-expss::vlookup(data1$Crosses,dict=data2,lookup_column = "Crosses",result_column = "BLUE")
  corBLUEs<-cor(data1$BLUE,data1$BLUE2,use="complete")
  return(corBLUEs)
}

cors(data1=BLUE3.9,data2=BLUE3.8) #0.97
cors(data1=BLUE3.9,data2=BLUE3.7) #0.98
cors(data1=BLUE3.9,data2=BLUE3.6) #0.94
cors(data1=BLUE3.9,data2=BLUE3.5) #0.82
cors(data1=BLUE3.9,data2=BLUE3)   #0.81

BLUE3.10$BLUE<-BLUE3.10$lsmean
cors(data1=BLUE3.9,data2=BLUE3.10)  #0.96

cors(BLUE3.8,BLUE3.5) #0.85
cors(BLUE3.8,BLUE3.6) #0.96
cors(BLUE3.8,BLUE3.7)  # 0.99
cors(BLUE3.8,BLUE3)  #0.85


## For the SAME deBLUP model as sommer, where popChk fixed; line+block+line:block+Crosses were random:
## asreml and sommer cors 
### for Yr2019 traits: 0.9611787 0.8991654 0.9814806 0.9935851 0.9743819 0.8527875
### for Yr2020 traits: 0.9339941 0.8484006 0.9863648 0.9932740 0.9924058 0.7861134

### Getting BLUEs, Across Years
# data$Trait<-data[,traits[i]]
# data<-data[order(data$Year),]
# mod7<-asreml(Trait~ line+block+popChk+Crosses,
#             random= ~Year:Crosses,
#             residual= ~ idv(units),
#             data = data, trace=TRUE,
#             maxiter=200,
#             na.action = na.method(y = "include",x="include"))

### Getting BLUEs ---- A different model
# mod6<-asreml(Trait~ popChk+Crosses,
#              random= ~line+block,
#              residual= ~ idv(units),
#              data = dataYr, trace=TRUE,
#              maxiter=200,
#              na.action = na.method(y = "include",x="include")) ## Higher AIC/BIC than mod3

### Getting d-BLUPs ---- A different model
# mod5<-asreml(Trait~ popChk,
#             random=~line+block+Crosses,
#             residual= ~ idv(units),
#             data = dataYr, trace=TRUE,
#             maxiter=200,
#             na.action = na.method(y = "include",x="include")) # Higher AIC/BIC than mod4

