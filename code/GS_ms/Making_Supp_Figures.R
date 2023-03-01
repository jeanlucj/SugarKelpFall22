#### Making Supp Figures

########### FIGURE 3 Plotting BVs correlation, and BVs by generations
WD<-"/local/workdir/mh865/GCA_SCA/" ## where all the BVs are stored
WD <- "/Users/jj332/Documents/GitRepo/SugarKelpBreedingMao/SugarKelpBreeding/TraitAnalyses201003/GS_ms/output_in_GS_ms/"

## To get the Average of Yr2019 + Yr2020 BLUPs
PlotTraits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades")
morTraits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

Averaging_BLUPs<-function(data=SingleBV,traits=PlotTraits){
  SingleBV$Year<-c(rep(2019,nrow(SingleBV)/2),rep(2020,nrow(SingleBV)/2))

  AllBLUPs<-NULL
  for (j in 1:length(traits)){
    SingleBV$trt<-SingleBV[,traits[j]]
    PlotBLUPMean<-aggregate(SingleBV$trt,list(SingleBV$Indiv),FUN=mean)  ### Just used mean BLUPs from Yr2019 and Yr2020!!!!
    rownames(PlotBLUPMean)<-PlotBLUPMean$Group.1
    PlotBLUPMean$trt<-rep(traits[j],nrow(PlotBLUPMean))
    AllBLUPs<-rbind(AllBLUPs,PlotBLUPMean)
  }
  ### Long to wide
    BLUPdfwide<-tidyr::spread(AllBLUPs,trt,x)
  return(BLUPdfwide)
}

### Plot level  A866  BVs
trtlevel<-"PlotLevel"
load(paste0(WD,"BVs_UsingYr1920_h2_SE_GenetCor_GScorPredA866Indiv_",trtlevel,".Rdata"))
plotBLUPs<-Averaging_BLUPs(data=SingleBV,traits=PlotTraits)
  rm(SingleBV)

### load Yr21 SP BVs
load(paste0(WD,"BVs_UsingYr1920_h2_SE_GenetCor_GScorPredA866Indiv_",trtlevel,".Rdata"))
Yr21BVplot<-Yr21BVfromComb
  rm(trtlevel)
  rm(Yr21BVfromComb)

### Indi level A866 BVs
trtlevel<-"IndiLevel"
load(paste0(WD,"BVs_UsingYr1920_h2_SE_GenetCor_GScorPredA866Indiv_",trtlevel,".Rdata"))
IndiBLUPs<-Averaging_BLUPs(data=SingleBV,traits=morTraits)
  rm(SingleBV)

### load Yr21 SP BVs
load(paste0(WD,"BVs_UsingYr1920_h2_SE_GenetCor_GScorPredA866Indiv_",trtlevel,".Rdata"))
Yr21BVIndiv<-Yr21BVfromComb

### Merge Plot and Indi level BVs for the 866 Crosses
allBLUPs<-merge(plotBLUPs,IndiBLUPs,by.x="Group.1",by.y="Group.1",all.x=TRUE)
rownames(allBLUPs)<-allBLUPs$Group.1
allBLUPs<-allBLUPs[,!colnames(allBLUPs)=="Group.1"]

### Yr21 SP BVs
Yr21BVall<-merge(Yr21BVplot,Yr21BVIndiv,by.x="row.names",by.y="row.names",all.x=TRUE)
rownames(Yr21BVall)<-Yr21BVall$Row.names
Yr21BVall<-Yr21BVall[,!colnames(Yr21BVall)=="Row.names"]

#### change the traits name
traits<-c(PlotTraits,morTraits)
traitsorder<-c("WWP","PDW","DWpM","AshOnly","AshFDwPM","BDns","BLen","BMax","BThk","SLen","SDia") ### Change it here!!!!!!!!!
for (i in 1:length(traits)){
  colnames(allBLUPs)<-ifelse(colnames(allBLUPs)==traits[i],traitsorder[i],colnames(allBLUPs))
  colnames(Yr21BVall)<-ifelse(colnames(Yr21BVall)==traits[i],traitsorder[i],colnames(Yr21BVall))
}

Yr21BVall<-Yr21BVall[,colnames(allBLUPs)]

#### This is for later to pull out the Crosses from each Year form the raw data
datafdr2<-paste0(WD,"OneTime192021/data/")
datafdr2 <- "/Users/jj332/Documents/GitRepo/SugarKelpBreedingMao/SugarKelpBreeding/TraitAnalyses201003/GS_ms/data_in_GS_ms/"
load(paste0(datafdr2,"3yrs_data_plotformated_indivNot_12012021.Rdata"))    ## Crosses in Yr21; #load: dataNHpi3yrs_C_Ash

#### Compare that to the A866 Individuals. RM the 4 checks
load(paste0(datafdr2,"outCovComb4_Mix_Conden_0527_2021_866Indiv.Rdata"))  #A866 load: outCovComb4_MixOrder

rownames(allBLUPs)[!rownames(allBLUPs)%in%rownames(outCovComb4_MixOrder)] ## Four checks
allBLUPs<-allBLUPs[rownames(allBLUPs)%in%rownames(outCovComb4_MixOrder),] ## Only kept those 866 in the A_866 list

allBLUPs<-allBLUPs[match(rownames(outCovComb4_MixOrder),rownames(allBLUPs)),]

#### Defining different categories for 866 Indiv (not including Yr21 SPs yet!!)
  Fndrs<-rownames(allBLUPs)[!(grepl("FG",rownames(allBLUPs))|grepl("MG",rownames(allBLUPs)))]
  FGMG_Fndr<-rownames(allBLUPs)[!(grepl("UCONN",rownames(allBLUPs)) | grepl("x",rownames(allBLUPs))) & (grepl("FG",rownames(allBLUPs)) | grepl("MG",rownames(allBLUPs)))]
  FGMG_UCONN<-rownames(allBLUPs)[grepl("UCONN",rownames(allBLUPs)) & !grepl("x",rownames(allBLUPs))]
  AllCrosses<-rownames(allBLUPs)[grepl("x",rownames(allBLUPs))]

  Crosses2019<-unique(as.character(dataNHpi3yrs_C_Ash[dataNHpi3yrs_C_Ash$Year==2019,]$Crosses))
  Crosses2020<-unique(as.character(dataNHpi3yrs_C_Ash[dataNHpi3yrs_C_Ash$Year==2020,]$Crosses))
  Crosses2021<-unique(as.character(dataNHpi3yrs_C_Ash[dataNHpi3yrs_C_Ash$Year==2021,]$Crosses))
  Crosses1920uniq<-unique(c(Crosses2019,Crosses2020))
  Crosses1920Common<-c(Crosses2019,Crosses2020)[duplicated(c(Crosses2019,Crosses2020))]
  Crosses2021Common<-c(Crosses2020,Crosses2021)[duplicated(c(Crosses2020,Crosses2021))]

  # CrossesSelf<-dataNHpi$Crosses[dataNHpi$isSelf==1]

  allBLUPs$Category<-NA

  allBLUPs$Category[rownames(allBLUPs)%in%Crosses2019 & !rownames(allBLUPs)%in%Crosses1920Common]<-"SP_2019"
  allBLUPs$Category[rownames(allBLUPs)%in%Crosses1920Common]<-"SP_2019"

  allBLUPs$Category[rownames(allBLUPs)%in%Crosses2020 & !rownames(allBLUPs)%in%Crosses1920Common & !rownames(allBLUPs)%in%Crosses2021Common]<-"SP_2020"

  allBLUPs$Category[rownames(allBLUPs)%in%Crosses2021 & !rownames(allBLUPs)%in%Crosses2021Common]<-"SP_2021"
  allBLUPs$Category[rownames(allBLUPs)%in%Crosses2021Common]<-"SP_2021"

  allBLUPs$Category[rownames(allBLUPs)%in%Fndrs]<-"Founders"
  allBLUPs$Category[rownames(allBLUPs)%in%FGMG_Fndr]<-"GP_from_Founders"
  allBLUPs$Category[rownames(allBLUPs)%in%FGMG_UCONN]<-"GP_from_Farm"

  #Fndrs, 104; FGMG_Fndr, 439; FGMG_UCONN, 78; AllCrosses, 329
  hasDesc <- sapply(Fndrs, function(n) grep(n, rownames(allBLUPs)[rownames(allBLUPs)%in%AllCrosses])) # grep the fndr pattern in progenySP string
  nDesc <- sapply(hasDesc, length)  #sapply: function applied to each element of a vector
  nDesc
  unique(allBLUPs$Category)

  allBLUPs$Category[which(nDesc==0)]<-"Founder_SP_without_progeny"

  ### RM the allBLUPs Category that are NAs  # There are 13 of crosses uncategorized
  ### FROM JLJ: I had 1 cross with Category NA: SA18-CB-5-FG1xSA18-CB-10-MG1
  allBLUPs<-allBLUPs[!is.na(allBLUPs$Category),]

  ############ I am merging the Yr2021 SPs combined BVs to the A866 allBLUPs file
  Yr21BVall$Category<-"SP_2021"
  allBLUPs0<-rbind(allBLUPs,Yr21BVall)

  ###### Adding groups to the allBLUPs0
  UCONN_GPto_SP<-rownames(allBLUPs0)[grepl("UCONN",rownames(allBLUPs0)) & grepl("x",rownames(allBLUPs0))]

  allBLUPs0$Group<-"SPs"
  allBLUPs0$Group[rownames(allBLUPs0)%in%UCONN_GPto_SP]<-"SPs_used_GP_from_Farm"  ###  This will be true if added Yr2021 SPs
  allBLUPs0$Group[allBLUPs0$Category%in%c("GP_from_Founders","GP_from_Farm")]<-"GPs"

  ####################################################################
  ###  JLJ to make boxplots
  library(ggplot2)
  p <- ggplot(allBLUPs0, aes(x=Category, y=DWpM)) +
    geom_boxplot(notch=T)
  p
  # This looks ugly because the category names are so long
  # Change the names
  allBLUPs0$shrtCat <- NA
  shrtCat <- c("FndSP+", "FndSp-", "GP_Fnd", "SP19_Fnd", "SP20_Fnd", "SP21", "GP_Frm")
  names(shrtCat) <- unique(allBLUPs0$Category)
  allBLUPs0$shrtCat <- shrtCat[allBLUPs0$Category]
  ggplot(allBLUPs0, aes(x=shrtCat, y=DWpM)) +
    geom_boxplot(notch=T) + scale_x_discrete(limits=c("FndSp-", "FndSP+", "GP_Fnd", "GP_Frm", "SP19_Fnd", "SP20_Fnd", "SP21"))

  allBLUPs0$shrtCat2 <- allBLUPs0$shrtCat
  allBLUPs0$shrtCat2[allBLUPs0$Group == "SPs_used_GP_from_Farm"] <- "SP21_Frm"
  allBLUPs0$shrtCat2[allBLUPs0$shrtCat2 == "SP21"] <- "SP21_Fnd"

  pdf("Supp Fig 3 boxplots.pdf")
  ggplot(allBLUPs0, aes(x=shrtCat2, y=DWpM)) +
    geom_boxplot(notch=T) +
    scale_x_discrete(limits=c("FndSp-", "FndSP+", "GP_Fnd", "GP_Frm", "SP19_Fnd", "SP20_Fnd", "SP21_Fnd", "SP21_Frm")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16)) +
    theme(axis.title.y = element_text(size=16)) +
    xlab("") + ylab("DWpM Breeding Values")
  dev.off()

#### Supp Fig 2 !!!
  data<-as.data.frame(allBLUPs0)
  data$Trait<-data$DWpM

  yr<-"Yr1920_AddYr21BV"
  diphap<-"MixedPloidy"
  traitname<-"DWpM"
  library(ggplot2)
  pdf(paste0(datafdr,"Supp Fig 2",traitname,"vs Individual",yr,"_PhotoScore23_with_",diphap,"_Supp_Fig_2_012022.pdf"))

  ggplot(data=data,mapping=aes(x=1:nrow(data),y=Trait,color=Category,shape=Group))+
    geom_point()+
    scale_color_manual(name="Category",
                       values=c("Founders"="black",
                                "Founder_SP_without_progeny"="grey",
                                "GP_from_Founders"="darkgoldenrod2",
                                "SP_2019"="blue1",

                                "SP_2020"="purple",

                                "SP_2021"="brown",

                                "GP_from_Farm"="tomato1"))+
    scale_shape_manual(values=c("SPs"=15,"SPs_used_GP_from_Farm"=3,"GPs"=17))+
    labs(x="Individual",y=paste0("BVs for ",traitname))+
    theme_bw()
  dev.off()

#### Supp Fig 1 !!!!!  Cor plot of SPs from Yr19 and Yr20 (n=232 plots)
SP_allBLUPs<-allBLUPs[rownames(allBLUPs)%in%AllCrosses,!colnames(allBLUPs)%in%c("Category","Group")]
pdf(paste0(datafdr,"Supp Fig 1_CorrPlot_SPplots",yr,"_PhotoScore23_With_",diphap,"_0120_2021.pdf"))
corrplot::corrplot.mixed(cor(SP_allBLUPs,use="complete"), diag="n", tl.cex=0.6, tl.col=1)
dev.off()


#### Supp Fig 3 !!!!!!  Still needs to find the newer script!!!!
SP_allBLUPs<-SP_allBLUPs[order(-SP_allBLUPs$DWpM),]
  dim(SP_allBLUPs)

plotBLUPs<-SP_allBLUPs
  head(plotBLUPs)
library(ggplot2)

### Had GPs biomass, being used in making Crosses this year

Plot_Cross_Link<-read.csv(paste0(datafdr,"Plot_num_Cross_Link_02052021.csv"),sep=",",header=T)

library(expss)
plotBLUPs$Isolate_Sorus<-vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Isolate_Sorus_2020_Need_Update",lookup_column = "Crosses")

plotBLUPs$Crossed2020<-vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Used_in2020Cross",lookup_column ="Crosses" )

## AshFDwPM

pdf(paste0(datafdr,"Supp Fig 3_Histogram_of_SP_BVs_which_were_selected.pdf"))
data<-plotBLUPs

data$Trait<-data$DWpM ## !!!
histo<-ggplot(data,aes(x=Trait))+
  geom_histogram(binwidth=0.01)

# segment_data = data.frame(
#   x = c(data$Trait[data$Isolate_Sorus=="Yes"]),
#   xend = c(data$Trait[data$Isolate_Sorus=="Yes"]),
#   y = rep(0,length(data$Trait[data$Isolate_Sorus=="Yes"])),
#   yend = rep(5, length(data$Trait[data$Isolate_Sorus=="Yes"]))
# )  ### red, the Plots were being isolated sorus tissue

segment_data2 = data.frame(
  x = c(data$Trait[data$Crossed2020=="y"]),
  xend = c(data$Trait[data$Crossed2020=="y"]),
  y = rep(0,length(data$Trait[data$Crossed2020=="y"])),
  yend = rep(4, length(data$Trait[data$Crossed2020=="y"]))
)  ### green, the 12 plots had GPs used for making crosses in Yr2020 to make Yr21 SPs. Some of these Yr21 SPs did not have useful data (as seen from allBLUPs0$Group=="SPs_used_GP_from_Farm")

plot<-histo +

  geom_segment(data=segment_data2,aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="green")+
  labs(x="DwPM",y="Plot Count")

#  geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+

print(plot)
dev.off()
