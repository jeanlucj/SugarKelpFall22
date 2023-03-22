library(tidyverse)
#Import data
dartagMrk <- read_csv("/Users/yaoguangli/Downloads/SugarKelpFall22/data/Report_DSacc22-7080_SNP.csv",
                      skip = 6,na = "-")

# Some curating of column names, copied from QualCtrlDArTag.R
cnames <- colnames(dartagMrk)
cnames <- gsub("Old", "OLD", cnames, fixed=T)
cnames <- gsub("G-", "G", cnames, fixed=T)
cnames <- gsub("UCONN-S", "UCONN-", cnames, fixed=T)
cnames <- gsub("UCONN-", "UCONN-S", cnames, fixed=T)
cnames <- gsub("SFG", "FG", cnames, fixed=T) #Added by Yaoguang
cnames[cnames == "SA-CB-7-MG3"] <- "SA18-CB-7-MG3"
missingS <- substr(cnames, 8, 8) == "-"
addS <- paste0(substr(cnames[missingS], 1, 8), "S", substr(cnames[missingS], 9, 100))
cnames[missingS] <- addS
colnames(dartagMrk) <- cnames

###Solve the duplicated sample names

duplicated_names <- duplicated(colnames(dartagMrk), fromLast = FALSE)
dartagMrk_DupCName <- dartagMrk %>%
  rename("SL20-MB-S2-FG3_2" = 179, "SL19-UCONN-S90-FG1_2" = 188)
colnames(dartagMrk_DupCName)

###change 2 to 0.5 for the marker data

dartagMrk_0.5Het <- dartagMrk_DupCName %>%
  mutate_at(c(17:204),funs(replace(., . > 1, 0.5)))

dartagMrk_0.5Het$`SL18-OI-S15-MG3` #randomly check one sample's data


##############Yanguang tries to make G matrix using AGHmatrix-----------start---------------
# Calculate relationship matrix G matrix using package AGHmatrix
# AGHmatrix Tutorial https://cran.r-project.org/web/packages/AGHmatrix/vignettes/Tutorial_AGHmatrix.html
# install.packages("devtools")
# devtools::install_github("rramadeu/AGHmatrix")
library(AGHmatrix)
MrK_OrigData <- dartagMrk_DupCName
# check the samples (columns) with >80% NAs
MrK_RM80NA <- MrK_OrigData[, which(colMeans(!is.na(MrK_OrigData)) > 0.8)]
#Remove SL20-MB-S2-FG3-2,SL20-MB-S6-MGOLD12, SL20-MB-S11-FGOLD1
MrK_RMSample <- select(MrK_RM80NA, -c("SL20-MB-S2-FG3_2","SL20-MB-S11-FGOLD1"))
#Make the matrix
MrK_Mat <- select(MrK_RMSample, -c(2:16))
# 0 = Hom Ref, 1 = Hom Alt, 2 = Het, change to 0 = Hom Ref, 2 = Hom Alt, 1 = Het
MrK_Mat[MrK_Mat==1] <- 3
MrK_Mat[MrK_Mat==2] <- 1
MrK_Mat[MrK_Mat==3] <- 2
#Transpose the tibble
MrK_Trans <- MrK_Mat
MrK_ready<- MrK_Trans %>%
  gather(var, val, 2:ncol(MrK_Trans)) %>%
  spread(MarkerName, val)

#Covert tibble to matrix
row_name = MrK_ready$var
MrK_ready %<>% select(-var) %>% as.matrix
rownames(MrK_ready) = row_name

MrK <- MrK_ready
str(MrK)


#Computing the additive relationship matrix based on VanRaden 2008
#VanRaden, P. 2008 Efficient methods to compute genomic predictions. Journal of dairy science 91, 4414–4423.
G_VanRadenPine <- Gmatrix(SNPmatrix=MrK, missingValue=-9,
                          maf=0.05, method="VanRaden")

#Computing the dominance relationship matrix based on Vitezica 2013
#Vitezica, ZG, et al. 2013 On the additive and dominant variance and covariance of individuals within the genomic selection scope. Genetics 195, 1223–1230
G_VitezicaPine <- Gmatrix(SNPmatrix=MrK, missingValue=-9,
                          maf=0.05, method="Vitezica")

##############Yanguang tries to make G matrix using AGHmatrix-----------End---------------

##########Yaoguang 3/15/2023 updates----------------------------------------
####JL function of calculating the relationship matrix, debugging----------
### Calculate the relationship matrix for GPs
calcGenomicRelationshipMatrix <- function(locusMat, ploidy=2){
  if (!any(ploidy == 1:2)) stop("Ploidy must be 1 or 2")
  freq <- colMeans(locusMat) / ploidy
  locusMat <- scale(locusMat, center=T, scale=F)
  return(tcrossprod(locusMat) / sum(ploidy*freq*(1-freq)))
}

tst <- dartagMrk_DupCName %>%
  mutate_at(c(17:204),funs(replace(., . > 1, 0.5)))

tst_Locus <- as.matrix(tst[,17:204])
tst_ind <- tst_Locus[ ,colMeans(is.na(tst_Locus))<0.2]
tst_marker <- tst_ind[rowMeans(is.na(tst_ind))<0.001, ]

tst_correlat <- calcGenomicRelationshipMatrix(t(as.matrix(tst_marker)), ploidy=1)
tst_correlat
######because of too many NAs in the matrix, the tcrossprod returns all NAs in the relationship matrix
####Need to clean the marker data file before calculate the relationship matrix



