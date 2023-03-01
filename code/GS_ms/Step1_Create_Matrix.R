# Estimate heritabilities with more power (w/o response to E)
# Pedigree tracking for measures of relatedness
# Estimate line and block affect using checks (C-1 and C-2 etc.)
# Estimated breeding value 
# Use called SNPs for marker-based relationship matrix
rm(list=ls())
library(here)
here()
#"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding"

dataDir<-here("TraitAnalyses201003/data")

  ls()
library(magrittr)
library(rrBLUP)

### Manually !!! Make 2020 Yr plotNo unique from 2019 Yr plotNo
###!!!!!!!!!!!!!!!!!!!!

#### I. Formatting Yr2019 and Yr2020 data
dataAll <-read.csv(paste0(dataDir,"/Allplotdata_03.08.2021.csv"),sep=",",header=TRUE) ##!!!
  str(dataAll)  
  
### Adding data with Ash ---- 0309_2021
CombineTrait<-function(data=dataAll,TraitName1="Ash..Celignis.",TraitName2="Ash..WHOI."){
  trait<-c()
  for (i in 1:nrow(data)){
    trait1<-data[,TraitName1][i]
    trait2<-data[,TraitName2][i]
    if (is.na(trait1)){
      trait[i]<-trait2
    }else if(!is.na(trait1) & !is.na(trait2)){
      trait[i]<-mean(c(trait1,trait2))
    }else if (!is.na(trait1) & is.na(trait2)){
      trait[i]<-trait1
    }
  }
  return(trait)
}

dataAll$Ash<-CombineTrait(data=dataAll,TraitName1 = "Ash..Celignis.",TraitName2 = "Ash..WHOI.")  
dataAll$C_N1<-dataAll$Carbon..Celignis./dataAll$Nitrogen..Celignis. 

dataAll$C_N_Combine<-CombineTrait(data=dataAll,TraitName1="C_N1",TraitName2="C.N..UCONN.")

dataAll$AshFDwPM<-(dataAll$wetWgtPerM*dataAll$percDryWgt/100)*(1-(dataAll$Ash/100))
dataAll$popChk <- ifelse(substr(dataAll$plotNo, 1, 1) == "Z", substr(dataAll$plotNo, 1, 2), "ES")  # Checks VS ES
  ls()
  head(dataAll)
  dataAll[1:10,] 
  dim(dataAll)  #616
########################### A. Plot level data input, estimate pedNH
dataNHpi<-dataAll[dataAll$Region=="GOM",] ## !!! RM SNE, 530 rows
  dim(dataNHpi)
dataNHpi <- dataNHpi[order(dataNHpi$plotNo),]  #### Plots in alphabetic order
  dim(dataNHpi)
  str(dataNHpi)    
dataNHpi<-dataNHpi[!dataNHpi$crossID=="Buffer",]  ## !!! RMed Buffer lines

#dataNHpi<-dataNHpi[!dataNHpi$PhotoScore==0,] ## !!! RM PhotoScore=0,  447 rows
dataNHpi<-dataNHpi[dataNHpi$PhotoScore>1,] ## !!! RM PhotoScore<=1
  dim(dataNHpi)   # 283
  colnames(dataNHpi)
dataNHpi19<-dataNHpi[dataNHpi$Year==2019,] ## 2019 data   !!!!!!
dataNHpi19_C<-dataNHpi19
dataNHpi19_C<-dataNHpi19_C[order(dataNHpi19_C$plotNo),]  ## Order plotNo alphabetically

dataNHpi20<-dataNHpi[dataNHpi$Year==2020,]   #
dataNHpi20_C<-dataNHpi20
dataNHpi20_C<-dataNHpi20_C[order(dataNHpi20_C$plotNo),]  ## Order plotNo alphabetically

dataNHpiBoth_C<-dataNHpi  ## Order plotNo alphabetically
dataNHpiBoth_C<-dataNHpiBoth_C[order(dataNHpiBoth_C$plotNo),] ## Order plotNo alphabetically

#save(dataNHpi19_C,dataNHpi20_C,dataNHpiBoth_C,file="dataNHpi_withChk_3_sets_PhotoScore23.rdata")
save(dataNHpi19_C,dataNHpi20_C,dataNHpiBoth_C,file=paste0(dataDir,"/dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))

########################### II. Run it Only Once !!!! Fndrs Amat
##############3 RM "checks"
dataNHpi_rmChk<-dataNHpi[!dataNHpi$crossID=="Check",]    # 250
FMGPs<-unique(c(as.character(dataNHpi_rmChk$femaPar),as.character(dataNHpi_rmChk$malePar)))
SPs<-strsplit(as.character(FMGPs), split="-", fixed=T)
CrossedSP<-unique(sapply(SPs, function(vec) paste(vec[1:3], collapse="-"))) #70 of the unique fndrs used in making crosses for GOM, adding SNE it is 90 of them
  str(FMGPs)  # 227 unique GPs being used in crosses
  head(CrossedSP)
  str(CrossedSP)
  tail(CrossedSP)

######### II. Estimaating fndersAmat
load("FarmCPU_GAPIT.Rdata")
  #write.csv(geno2$taxa,"Fndr_Crossed_genotyped_Dart1.csv")
  ls()
geno2<-geno[,-1]
rownames(geno2)<-geno$taxa
  geno2[1:5,1:6]
  dim(geno2)
rownames(geno2) <- paste0(substring(rownames(geno2), first=1, last=2), "18", substring(rownames(geno2), first=3))
geno2[geno2==0]=-1
geno2[geno2==1]=0
geno2[geno2==2]=1
  
geno2<-geno2[rownames(geno2)%in%CrossedSP,]  # 70 fnders made crosses, only 47 were genotyped
  #write.csv(rownames(geno2),"Fndrs_genotyped_Samples.csv")
fndrMrkData<-geno2
  
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T,shrink=TRUE) ## Add shrink per Deniz
fndrMrkDataImp <- mrkRelMat$imputed
mrkRelMat <- mrkRelMat$A  ####### !!! 1. Fnder Amat
  dim(mrkRelMat)  #47 x 47  !!!!!!!
  #write.csv(rownames(mrkRelMat),"Fndrs_used_in_Amat.csv")
  save(mrkRelMat,file="fndr_mrkRelMat.Rdata") 
###########################  
  
## RM checks to make the pedigree-biphasic matrix and CC-matrix
dataNHpi19<-dataNHpi19[!dataNHpi19$crossID=="Check",]  
dataNHpi20<-dataNHpi20[!dataNHpi20$crossID=="Check",]  
dataNHpiBoth<-dataNHpi[!dataNHpi$crossID=="Check",]
  dim(dataNHpi19) #122
  dim(dataNHpi20)  #128
  dim(dataNHpiBoth)  #250
  
####Steps in the below ## Part
###*** Pedigree matrix calculation was done in the "Redo_Pedigree.R" file
  
# 2.Calculate the CC matrix: relationship matrix
source("calcCCmatrixBiphasic.R")
biphasicPedNH<-read.csv("Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv.csv",sep=",",header=TRUE,row.names=1)

biphasicCCmat <- calcCCmatrixBiphasic(biphasicPedNH)  #### !!!!!!! Update into calcCCmatrixHaploid() ?????
rownames(biphasicCCmat) <- colnames(biphasicCCmat) <- rownames(biphasicPedNH)


#   HaploidCCcal <- calcCCmatrixHaploid(biphasicPedNH)  
#   biphasicCCmat<-HaploidCCcal$mccMat ### !!!!!!!
#   biphasichccMat<-HaploidCCcal$hccMat
###### This is to work the  the haploid level of matrices
#### Take out the progeny SPs from pedigree matrix, put it back in later
#### Founders marker data, make haploid level, one sample 2 rows: 0->0,0 1->0,1 2->1,1
#### Cal relationship matrix for the fnders
#### ----> Use the paper method to combine this Founder relationship gmat, and the GP sequenced amat and the PedNH


# 3. Combine the marker- with the pedigree- relationship matrix to make H matrix #### !!!!! How to do that?????
aMat <- 2 * biphasicCCmat
  dim(aMat)  ###### !!! 3. Pedigree based CC matrix

save(aMat,file="Re_order_aMat.Rdata")
  
# fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))  ### Founders, which() only keeps the TRUE ones
# gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))   ### GPs, the col 3 is NA
# spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
# nSp <- length(spRows)
#   length(fndRows) # 73; 78; 79
#   length(gpRows) #99; 156; 227
#   length(spRows) #122; 127; 244+1
#   length(fndRows)+length(gpRows)+length(spRows) #294; 361; 550

hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))
save(mrkRelMat,aMat,hMat,biphasicPedNH,biphasicCCmat,fndrMrkData, file=paste0("hMat_PedNH_CCmat_fndrMrkData_",yr,"_PhotoScore23_WithSGP_866_Reorder_Pedigree.rdata"))  ###### !!!!!!!



############# B. Individual datafile formating
### Manually correct the plotNo to be unique in 2019 and 2020
dataNHim <- read.csv("Indi_2019_2020_Compile_07202020_Edit_CorrectplotNo.csv",sep=",",header=TRUE)
  dim(dataNHim)
  str(dataNHim)
  tail(dataNHim)

dataNHim<-dataNHim[dataNHim$Region=="GOM",] ## !!! RM SNE, 530 rows
  dim(dataNHim)  #5046

dataNHim<-dataNHim[!dataNHim$crossID=="Buffer",]
  dim(dataNHim) #4996

dataNHim19_C<-dataNHim[dataNHim$Year==2019,] 
dataNHim19_C<-dataNHim19_C[order(dataNHim19_C$plotNo),]  ### Plots in Alphabetic order; Here the plot did not RM plots with photoScore < 2,3
  dim(dataNHim19_C)  # 2800 x 12
  tail(dataNHim19_C)

dataNHim20_C<-dataNHim[dataNHim$Year==2020,]
dataNHim20_C<-dataNHim20_C[order(dataNHim20_C$plotNo),]
  dim(dataNHim20_C)  # 2196 x 12
  tail(dataNHim20_C)
  
dataNHimboth_C<-dataNHim
dataNHimboth_C<-dataNHimboth_C[order(dataNHimboth_C$plotNo),]
  tail(dataNHimboth_C)
###dataNHim$plotNo <- gsub("S", "", dataNHim$plotNo, fixed=TRUE) ### The 2020 plotNo is "SXXX"; 2019 is just number
###dataNHim$plotNo[dataNHim$plotNo == ""] <- "C2-I"

save(dataNHim19_C,dataNHim20_C,dataNHimboth_C,file="dataNHim_withChk_3_sets_PhotoScore0123.rdata")
#### !!!! This individual dataset still HAS the photoscore < 2,3 plots


# #### *** Below is to make the biphasicPedNH file, which is now reordered and remade in Redo_Pedigree.R   
# ######################################################################  
# ### Make pedigree relationship matrix  !!!!!!!!!! RUN 3 TIMES for 3 data sets !!!!!! 
# 
# #1. Create the pedigree
# source("makeBiphasicPed.R")
# 
# # dataNHpi<-dataNHpi19  ######### !!!!
# # yr<-"2019"
# # 
# # dataNHpi<-dataNHpi20  ######### !!!!
# # yr<-"2020"
# 
# dataNHpi<-dataNHpiBoth ########## !!!!!
# yr<-"Both"
# 
# dataNHpi$Crosses<-as.factor(as.character(dataNHpi$Crosses))
# dataNHpi$femaPar<-as.character(dataNHpi$femaPar)
# dataNHpi$malePar<-as.character(dataNHpi$malePar)
# 
# dataNHpi<-dataNHpi[order(dataNHpi$plotNo),]  ### Order plots alphabetically, so that order match that in making the Z matrix
# 
 kelpNameColumns <- dataNHpi[, c("femaPar", "malePar")]
 kelpNameColumns <- kelpNameColumns[kelpNameColumns[,1] != "",]   #### ? Is this removing the Checks?
# 
# ### Need to remove the femaPar x malePar that's the same cross
 kelpNameColumns$RM<-paste0(kelpNameColumns[,1],"x",kelpNameColumns[,2])
 kelpNameColumns<-kelpNameColumns[!duplicated(kelpNameColumns$RM),]
 kelpNameColumns<-kelpNameColumns[,-3]
# 
# ##### This is first step, using the list of FG and MG used in crossing, for generating the dataNHpi
# biphasicPedNH <- makeBiphasicPed(kelpNameColumns, rownames(mrkRelMat)) #### The rownames(mrkRelMat) is to put the fndrs in the front in the pedigree matrix
#   head(kelpNameColumns)
#   head(rownames(mrkRelMat))
#   dim(biphasicPedNH)
#   tail(biphasicPedNH)
# write.csv(biphasicPedNH,paste0("biphasicPedNH_",yr,"_PhotoScore23_NoSGPs.csv") ) 
# 
# 
# ############
# ##### Run this only ONCE !!!!!
# fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))  ### Founders, which() only keeps the TRUE ones
# gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))   ### GPs, the col 3 is NA
# spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
# nSp <- length(spRows)
#   nSp
#   str(fndRows)
# 
#   
# ##### V I. Add additional GPs genotyped, wfndr, woutfndr  (Step 2)
# ### Add in the list of GPs being genotyped but are not in the crosses from 3 plates 
# 
 GPsequencedAll<-read.csv("3_plates_sequenced_names.csv",sep=",",header=T)
 ### GPsequenced1: NO checks Yet (unique list)
 GPsequenced1<-GPsequencedAll[!GPsequencedAll$Chosen.Reason=="Checks",]$Name_InPheno
 GPsequenced1<-unique(as.character(GPsequenced1))
 
 # ### GPsequenced2: checks Only (uqniue)
 # GPsequenced2<-GPsequencedAll[GPsequencedAll$Chosen.Reason=="Checks",]$Name_InPheno  
 # GPsequenced2<-unique(as.character(GPsequenced2))
 #"SL18-SF-13-MG1"
 
 GPsequenced<-GPsequenced1  ######## !!!!!
# 
#   
# #####{{}} 
# GPsequenced_fndr<-strsplit(GPsequenced,split="-", fixed=T) 
# GPsequenced_fndr<-sapply(GPsequenced_fndr, function(vec) paste(vec[1:3], collapse="-")) ## This is their fndr
# 
# ## This is GP already has fndr
# GPsequenced_w_fndr<-GPsequenced[which(GPsequenced_fndr%in%names(fndRows))]   ## GP_fndrs names in the fndRows
# #GPsequenced_fndr[which(GPsequenced_fndr%in%names(fndRows))] ## This is the fndrs of those GPs
# GPsequenced_w_fndr_unique<-GPsequenced_w_fndr[which(!GPsequenced_w_fndr%in%names(gpRows))] 
# GPsequenced_ws_fndr_unique<-strsplit(GPsequenced_w_fndr_unique,split="-", fixed=T) 
# GPsequenced_ws_fndr_unique<-sapply(GPsequenced_ws_fndr_unique,function(vec) paste(vec[1:3], collapse="-"))
# 
# ## This is GP without (wout) fndr 
# GPsequenced_wout_fndr<-GPsequenced[which(!GPsequenced_fndr%in%names(fndRows))]   
# GPsequenced_wouts_fndr<-GPsequenced_fndr[which(!GPsequenced_fndr%in%names(fndRows))] ## This is the fndrs of those GPs
# 
#   head(GPsequenced_wout_fndr)
#   head(GPsequenced_wouts_fndr)
#   tail(GPsequenced_wout_fndr)
#   tail(GPsequenced_wouts_fndr)
#   head(GPsequenced_w_fndr)
#   head(GPsequenced_ws_fndr_unique)
#   tail(GPsequenced_w_fndr)
#   tail(GPsequenced_ws_fndr_unique)
#   tail(biphasicPedNH)
#   dim(biphasicPedNH)
# 
# biphasicPedNH0<-biphasicPedNH 
# rowStart<-nrow(biphasicPedNH)+1
# 
# ### the Pedigree matrix for the GPs already has fndr  
# ### At the same time, RM those GPs that were already in the gpRows list
# 
# GPseq_wfndr_Ped<-cbind(rowStart:(rowStart-1+length(GPsequenced_w_fndr_unique)),fndRows[GPsequenced_ws_fndr_unique],NA)
# rownames(GPseq_wfndr_Ped)<-GPsequenced_w_fndr_unique
#   head(GPseq_wfndr_Ped)
#   tail(GPseq_wfndr_Ped)
#   dim(GPseq_wfndr_Ped)  
#   
# biphasicPedNH<-rbind(biphasicPedNH,GPseq_wfndr_Ped)
# rowStart2<-nrow(biphasicPedNH)+1
#   rowStart2
#   550+75
# GPsequenced_wouts_fndr_uniq<-unique(GPsequenced_wouts_fndr)
# GPseq_newfndr<-c(rowStart2:((rowStart2)-1+length(GPsequenced_wouts_fndr_uniq)))
# names(GPseq_newfndr)<-GPsequenced_wouts_fndr_uniq
# 
# ### the Pedi matrix part for the new Fndrs
# GPseq_newfndr_Ped<-cbind(GPseq_newfndr,0,0)
# 
# biphasicPedNH<-rbind(biphasicPedNH,GPseq_newfndr_Ped)
# rowStart3<-nrow(biphasicPedNH)+1
#   rowStart3
# 
# GPsequenced_wout_fndr<-unique(GPsequenced_wout_fndr)    
# GPseq_woutfndr_Ped<-cbind(rowStart3:(rowStart3-1+length(GPsequenced_wout_fndr)),GPseq_newfndr[GPsequenced_wouts_fndr],NA)
# rownames(GPseq_woutfndr_Ped)<-GPsequenced_wout_fndr
#   
# ### the Pedi matrix part for the GPs added new Fndrs
# biphasicPedNH<-rbind(biphasicPedNH,GPseq_woutfndr_Ped)
# rowStart4<-nrow(biphasicPedNH)+1
#   rowStart4
# 
# biphasicPedNH<-rbind(biphasicPedNH0,GPseq_wfndr_Ped,GPseq_newfndr_Ped,GPseq_woutfndr_Ped)  
#   dim(biphasicPedNH)
# ######### {{}}
#   
#   
# #####{} V.II Add the S plots into the pedigree  (Step 3, these are the GPs came out of the 2019 S-plots)
 summary<-read.csv("Active_culture_biomass_assesment_for_Mao_100120.csv",sep=",",header=T,na.strings=c("","NA"))  
   dim(summary)  
   head(summary)
   str(summary)
 dataAll <-read.csv("dataNHpi_ManuallyAddGrowth_9.29.20_Ireland_Data_AshFDWpM.csv",sep=",",header=TRUE) ##!!!
 
 fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))  ### Founders, which() only keeps the TRUE ones
 gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))   ### GPs, the col 3 is NA
 spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
 nSp <- length(spRows)
   length(fndRows)
   length(gpRows)
   nSp
   
 summary$GametophyteID<-as.character(summary$GametophyteID)  
 SPGPs<-summary[grep("UCONN-S",summary$GametophyteID),] 
 #SPGPs<-summary[summary$Generation=="SP_GP"&!is.na(summary$X..Crosses)&!summary$X..Crosses==0,] # RM those w/out biomass
   nrow(SPGPs)==length(unique(SPGPs$GametophyteID)) ### This MUST be TRUE to make the list of S_GPs UNIQUE
 SPGPs<-droplevels(SPGPs)
   dim(SPGPs)
   dim(SPGPs_ls)
   head(SPGPs)
 # SPGPs_ls<-as.data.frame(unique(SPGPs$GametophyteID))
 # colnames(SPGPs_ls)<-"GametophyteID"
 #   str(SPGPs_ls)
 #   head(SPGPs_ls)  
# 
# SPGPs0<-SPGPs
# SPGPs0$SPCross<-vlookup(SPGPs0$Parent,dict=dataAll,result_column="Crosses",lookup_column = "crossID")
# write.csv(SPGPs0,"SPGPs0.csv")
# 
# SPGPs<-SPGPs_ls
# SPGPs$Parent<-vlookup(SPGPs$GametophyteID,dict=SPGPs0,result_column = "Parent",lookup_column="GametophyteID")
# 
# #SPGPs$Parent<-vlookup(SPGPs$GametophyteID,dict=SPGPs0,result_column = "Parent",lookup_column="GametophyteID")
#   head(SPGPs)
#   dim(SPGPs)
# 
# library(expss)  
# SPGPs$SPCross<-vlookup(SPGPs$Parent,dict=dataAll,result_column="Crosses",lookup_column = "crossID")
#   head(SPGPs)
#   tail(SPGPs)
#   str(SPGPs)
# SPGPs$GametophyteID<-as.character(SPGPs$GametophyteID)  
# ## IF there is any SP_GP that did not have a plotlevel cross in the biphasicPed (its fndr, Plotcross)
#   ## That would make it more complicated. # Those NA SP_plot needs: fndr,plot,SP_GP\
#     #### SPGPs has the list of SP_GPs to be included
# 
# # ## RM the row with all NAs
# ####
# SPGPs$Still<-ifelse(SPGPs$GametophyteID%in%rownames(biphasicPedNH),"no","still")
# SPGPs2<-SPGPs[SPGPs$Still=="still",]
#   str(SPGPs2)
# SPGPs3<-SPGPs
# SPGPs<-SPGPs2    
# ####
# 
# #SP_GP<-SPGPs[unique(SPGPs$GametophyteID),]
# #SP_GP<-droplevels(SPGPs[unique(SPGPs$GametophyteID),])
#   dim(SP_GP)
#   head(SP_GP)
#   tail(SP_GP)
# #SP_GP<-SP_GP[apply(SP_GP, 1, function(v) !all(is.na(v))),]
# 
# SP_GP<-SPGPs
# 
# rowStart5<-nrow(biphasicPedNH)+1
#   rowStart5
# SP_GP$SPCross<-as.character(SP_GP$SPCross) ## !!! Need to be character for the below to work
# SP_GP$withPlotPedi<-vlookup(SP_GP$SPCross,dict=biphasicPedNH,result_column = 1,lookup_column = "row.names") ## !! Look up plotPedi row number
# 
# SPGP_withPlot<-SP_GP[!is.na(SP_GP$withPlotPedi),]
# SPGP_withPlot<-droplevels(SPGP_withPlot)
#   dim(SPGP_withPlot)
# 
# ### 1.1  # number of elements: in SPGP_woithPlot
# SPGP_withPlotpedi<-cbind(rowStart5:(rowStart5-1+nrow(SPGP_withPlot)),spRows[as.character(SPGP_withPlot$SPCross)],NA)
#   tail(SPGP_withPlotpedi)
#   str(SPGP_withPlotpedi)
#   str(SPGP_withPlot$GametophyteID)
# rownames(SPGP_withPlotpedi)<-SPGP_withPlot$GametophyteID
# 
# rowStart6<-SPGP_withPlotpedi[nrow(SPGP_withPlotpedi),1]+1
#   rowStart6
#   
# SP_GP_NA<-SP_GP[is.na(SP_GP$withPlotPedi),]
# SP_GP_NA<-droplevels(SP_GP_NA)
#   SP_GP_NA
# SP_GP_NA_GPs<-str_split_fixed(string=as.character(SP_GP_NA$SPCross), "x", 2)
# SP_GP_NA$FG<-SP_GP_NA_GPs[,1]
# SP_GP_NA$MG<-SP_GP_NA_GPs[,2]
# 
# SP_GP_NA$withFGPedi<-vlookup(SP_GP_NA$FG,dict=biphasicPedNH,result_column =1, lookup_column = "row.names")
# SP_GP_NA$withMGPedi<-vlookup(SP_GP_NA$MG,dict=biphasicPedNH,result_column =1, lookup_column = "row.names")
#     SP_GP_NA
#     
# ### Manually checked: All has the GPs included in the biphasicPed
# SP_GP_NA_parent<-unique(SP_GP_NA[,colnames(SP_GP_NA)%in%c("SPCross","FG","MG")])
#     SP_GP_NA_parent
#     
# ###1.2 Only had 1 element
# SP_GP_NA_parent_pedi<-cbind(rowStart6,gpRows[as.character(SP_GP_NA_parent$FG)],gpRows[as.character(SP_GP_NA_parent$MG)])    
#   SP_GP_NA_parent
# rownames(SP_GP_NA_parent_pedi)<-SP_GP_NA_parent$SPCross
#   SP_GP_NA_parent_pedi
#   
# # SP_GP_NA_parent_pedi_v<-as.vector(SP_GP_NA_parent_pedi[,1])
# # names(SP_GP_NA_parent_pedi_v)<-rownames(SP_GP_NA_parent_pedi) # get the vector of parents SP plots names
# 
# ###1.3 number of elements in SP_GP_NA
# rowStart7<-SP_GP_NA_parent_pedi[nrow(SP_GP_NA_parent_pedi),1]+1
#   rowStart7
# SP_GP_NA_itself_pedi<-cbind(rowStart7:(rowStart7-1+nrow(SP_GP_NA)),SP_GP_NA_parent_pedi[rownames(SP_GP_NA_parent_pedi)==SP_GP_NA_parent$SPCross,1],NA)  # Only 1 extra
# rownames(SP_GP_NA_itself_pedi)<-SP_GP_NA$GametophyteID 
#   SP_GP_NA_itself_pedi
# 
# biphasicPedNH1<-biphasicPedNH
# 
# biphasicPedNH<-rbind(biphasicPedNH1,SPGP_withPlotpedi,SP_GP_NA_parent_pedi,SP_GP_NA_itself_pedi)
#   dim(biphasicPedNH)  #754
#   tail(biphasicPedNH)
# ##{}
# 
#   #biphasicPedNH1<-biphasicPedNH
#   #biphasicPedNH<-rbind(biphasicPedNH1,SPGP_withPlotpedi)
#   
#   ######## Adding GPs has biomass not in the biPhasic list
#   GPsequencedAll<-read.csv("AddGP_to_Pedi.csv",sep=",",header=T)
#   colnames(GPsequencedAll)<-"Name_InPheno"
#   head(GPsequencedAll)
#   dim(GPsequencedAll)
#   GPsequenced<-unique(as.character(GPsequencedAll$Name_InPheno))
#   dim(GPsequenced)
#   str(GPsequenced)
#   
#   # Run {{}} this part  
#   
#   
# write.csv(biphasicPedNH,"biphasiPedNH_addmoreGP_866.csv")  
#   
# # #### Checks Need to !!! run the {{}} part
# # Checks<-read.csv("Checks_2019_2020.csv",sep=",",header=T)
# #   head(Checks)
# # Checks$femaPar<-as.character(Checks$femaPar)  
# # Checks$malePar<-as.character(Checks$malePar)  
# # 
# # GPchecks<-unique(c(Checks$femaPar,Checks$malePar))  
# # GPsequenced<-GPchecks  ######## Use this to Run up to the {{}}
# # 
# # ####Used output, then Add the checks GPs wout Pedi
# # rowStart2??
# # GPchks_wFndr_Ped<-GPseq_wfndr_Ped
# # GPchks_noFndr<-GPsequenced_wout_fndr
# # 
# # GPchks_noFndr_Ped<-cbind(rowStart2:(rowStart2-1+length(GPchks_noFndr)),NA,NA) ##?
# # rownames(GPchks_noFndr_Ped)<-GPchks_noFndr
# # 
# # biphasicPedNH<-rbind(biphasicPedNH,GPchks_wFndr_Ped,GPchks_noFndr_Ped)
# # GPchks_AllGPs<-c(GPchks_wFndr_Ped[,1],GPchks_noFndr_Ped[,1])
# # 
# # #### Add all the checks SP plots in the Pedi
# # rowStart<-nrow(biphasicPedNH)+1
# # 
# # GPchks_SPlots<-cbind(rowStart:(rowStart-1+nrow(GPchks_AllGPs_Ped)), GPchks_AllGPs[Checks$femaPar], GPchks_AllGPs[Checks$malePar])
# # #### Rows not correct  !!!!!
# 
# ####*** End here: The above has been re-made in "Redo_Pedigree.R"
# ###############################################################
# ################

