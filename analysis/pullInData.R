library(tidyverse)

here::i_am("analysis/pullInData.R")

# Document what QualCtrl generates:
# badIndCall and badIndHet show which individuals have low call rate or high
# heterozygosity rate out of 188 individuals with DArTag
# indCallRate and indFreqHet have all the information: named vectors
# mrkCallRate and mrkFreqHet: named vectors
source(here::here("analysis", "QualCtrlDArTag.R"))

########################################
# Phenotypes
# This file came from downloading all phenotypes from trials where the
# name includes "FARM"
# Initially large matrix with individuals in rows and traits in columns
# Lots of missing values
phenoData <- read_csv(file=here::here("data", "AllFarmPhenotypes.csv"),
                      skip=3, guess_max=10000, na=c("", "NA", "N/A"))
phenoData <- phenoData %>% filter(observationLevel == "plot")
sumIsNA <- apply(phenoData, 2, function(v) sum(is.na(v)))
phenoData <- phenoData %>% select(which(sumIsNA<nrow(phenoData)))
phenoData <- phenoData %>% select(-contains("DbId"))
phenoData <- phenoData %>% select(-contains("program"))
# Remove individuals that have very little phenotypic data
# Total of 66 traits...
# "Presence of sorus tissue abs/pres|CO_360:0000253"
# "Successful release unsucc/succ|CO_360:0000320"
justTraits <- phenoData %>% select(contains("CO_360")) %>%
  select(-contains(c("253", "320")))
hasData <- apply(justTraits, 1, function(v) sum(!is.na(v))) > 2
# 975 plot-level observations
phenoData <- phenoData %>% dplyr::filter(hasData)
# Just the ones with OK photo scores
phenoData123 <- phenoData %>% dplyr::filter(`Photo score 0-3|CO_360:0000300` > 0)
# Just the ones with OK photo scores and from GOM
phenoData123G <- phenoData123 %>% dplyr::slice(-grep("SNE", studyName))

########################################
# Pedigree data
pedData <- read_tsv(file=here::here("data", "AllAccessionsPedigree.txt"),
                    na=c("", "NA", "N/A"))

########################################
# Manual curating
# SA18-CB-FG5 can't be right because it should have come from an SP so it should
# be something like SA18-CB-SX-FG5
pedData <- pedData %>% filter(Accession != "SA18-CB-FG5")
# These "SL19-UCONN-S156-FG3"  "SL19-UCONN-S158-FG3"  "SL19-UCONN-S23-FG4"
# Have exact matches in the database: not sure why they are not in the pedigree
# So we should add them to the pedigree
addToPedData <-  tibble(
  Accession=c("SL19-UCONN-S156-FG3", "SL19-UCONN-S158-FG3", "SL19-UCONN-S23-FG4"),
  Female_Parent=c("SL19-UCONN-S156", "SL19-UCONN-S158", "SL19-UCONN-S23"),
  Male_Parent=c("SL19-UCONN-S156", "SL19-UCONN-S158", "SL19-UCONN-S23"),
  Cross_Type=rep("self, 3"))
pedData <- dplyr::bind_rows(pedData, addToPedData)

# NOTE there are a lot of duplicate accessions in the pedData
pedData$Accession %>% duplicated %>% sum # 1123 (!)

# Let's look at these duplications: if the Accession is duplicated the
# remainder of the row is too.  So just remove duplicates
dupAcc <- pedData %>% dplyr::filter(Accession %>% duplicated) %>%
  pull(Accession) %>% unique
dupPD <- pedData %>% dplyr::filter(Accession %in% dupAcc)
disPD <- dupPD %>% (dplyr::distinct)

pedData <- pedData %>% dplyr::filter(!duplicated(Accession))

# There are some SPs that came from mixed GPs that have a -Mixed suffix
# This makes them look like they are GPs (3 dashes) change to Mixed (no dash)
hasDashMx <- grep("-Mixed", pedData$Accession, fixed=T)
fixMixed <- function(s){
  whereDash3 <- gregexpr("-", s, fixed=T)[[1]][3]
  return(paste0(substr(s, 1, whereDash3-1), "Mixed"))
}
pedData$Accession[hasDashMx] <- sapply(pedData$Accession[hasDashMx], fixMixed)

# There are some SPs whose pedigrees are given in terms of their grand-parental
# SPs.  I don't want that.
nDash <- pedData$Accession %>%
  sapply(function(s) gregexpr("-", s, fixed=T)[[1]] %>% length)
nDashFP <- pedData$Female_Parent %>%
  sapply(function(s) gregexpr("-", s, fixed=T)[[1]] %>% length)
nDashMP <- pedData$Male_Parent %>%
  sapply(function(s) gregexpr("-", s, fixed=T)[[1]] %>% length)
SPtoSP_F <- pedData %>% dplyr::filter(nDash==2 & nDashFP==2)
SPtoSP_M <- pedData %>% dplyr::filter(nDash==2 & nDashMP==2)
# This is a mess
# Eliminate ones where the parent starts with SL18-UCONN-G
# Eliminate ones where the parent starts with LIS
# Eliminate ones where the parent contains ME
# The same sets with male and female SP parents
SPtoSP_F <- SPtoSP_F %>%
  dplyr::slice(-grep("SL18-UCONN-G", Accession, fixed=T)) %>%
  dplyr::slice(-grep("SL18-UCONN-G", Female_Parent, fixed=T)) %>%
  dplyr::slice(-grep("LIS", Female_Parent, fixed=T))

SPtoSP_M <- SPtoSP_M %>%
  dplyr::slice(-grep("SL18-UCONN-G", Male_Parent, fixed=T)) %>%
  dplyr::slice(-grep("ME", Male_Parent, fixed=T)) %>%
  dplyr::slice(-grep("PI", Male_Parent, fixed=T))

# This works on pedData directly... So watch out...
for (femPar in SPtoSP_F$Female_Parent %>% unique){
  whRow <- which(SPtoSP_F$Female_Parent == femPar)
  cnt <- 1
  for (row in whRow){
    # Change nomenclature if it's a Mixed SP
    addMx <- ifelse(gregexpr("Mixed", SPtoSP_F$Accession[row])[[1]] %>% length > 0, "Mixed", "")
    pedDataRow <- which(pedData$Accession == SPtoSP_F$Accession[row])
    pedData <- pedData %>% slice(-pedDataRow)
    pedData <- pedData %>% dplyr::add_row(
      Accession=SPtoSP_F$Accession[row],
      Female_Parent=paste0(femPar, "-FG", addMx, cnt),
      Male_Parent=paste0(femPar, "-MG", addMx, cnt),
      Cross_Type="biparental",
      .before=pedDataRow
    )
    pedData <- pedData %>% dplyr::add_row(
      Accession=paste0(femPar, "-FG", addMx, cnt),
      Female_Parent=femPar,
      Male_Parent=femPar,
      Cross_Type="self",
      .before=pedDataRow
    )
    pedData <- pedData %>% dplyr::add_row(
      Accession=paste0(femPar, "-MG", addMx, cnt),
      Female_Parent=femPar,
      Male_Parent=femPar,
      Cross_Type="self",
      .before=pedDataRow
    )
    cnt <- cnt+1
  }
}

# Separate out SPs from GPs
nDash <- pedData$Accession %>%
  sapply(function(s) gregexpr("-", s, fixed=T)[[1]] %>% length)
pedDataSP <- pedData %>% filter(nDash == 2)
pedDataGP <- pedData %>% filter(nDash > 2)

# NOTE: I have to do some curation using the markers on the DArTag to eliminate
# the ones that look like they are mixtures

########################################
# Marker relationship matrices
# NOTE: these relationship matrices were downloaded directly from sugarkelpbase
# I need to specify what parameters were used in terms of missing and MAF
fndRelMat <- read_tsv(file=here::here("data", "DArTSeqRelationshipMatrix029910.tsv"), na=c("", "NA", "N/A"))
fndAccInPed <- fndRelMat$stock_uniquename %in% pedData$Accession
print(sum(fndAccInPed)) # All 179 founders in pedigree

tagRelMat <- read_tsv(file=here::here("data", "DArTagRelationshipMatrix029910.tsv"), na=c("", "NA", "N/A"))
tagAccInPed <- tagRelMat$stock_uniquename %in% pedData$Accession
print(sum(tagAccInPed)) # 175 out of 182 DArTag'ed accessions in pedigree
tagRelMat$stock_uniquename[!tagAccInPed]
# [1] "SL18-SF-S13-SFG4"     "SL18-SF-S13-SFG5"     "SL18-SF-S13-SFG6"
# [4] "SL19-UCONN-S90-FG1-2" "SL20-MB-S11-FGOLD1"   "SL20-MB-S2-FG3-2"
# [7] "SL20-MB-S6-MGOLD12"
# Remove the SFG names
tagRelMat$stock_uniquename <- gsub("SFG", "FG", tagRelMat$stock_uniquename)

# SL19-UCONN-S90-FG1-2. Both SL19-UCONN-S90-FG1 and SL19-UCONN-S90-FG2 are in
# the database but do not have exactly the same marker profile as SL19-UCONN-S90-FG1-2.  So we don't know what it is.
# SL20-MB-S2-FG3-2 is 100% similar to SL20-MB-S2-FG3.  Match.
# So drop SL20-MB-S2-FG3-2
# SL20-MB-S6-MGOLD12 is the same as SL20-MB-S6-MG12.
# So drop SL20-MB-S6-MGOLD12
# SL20-MB-S11-FGOLD1 is the same as SL20-MB-S11-FG1
# So drop SL20-MB-S11-FGOLD1
keepInTagData <- which(!(tagRelMat$stock_uniquename %in% c("SL19-UCONN-S90-FG1-2", "SL20-MB-S11-FGOLD1", "SL20-MB-S2-FG3-2", "SL20-MB-S6-MGOLD12")))
tagRelMat <- tagRelMat %>% slice(keepInTagData) %>%
  select(c(0, keepInTagData)+1)
tagAccInPed <- tagRelMat$stock_uniquename %in% pedData$Accession
print(sum(tagAccInPed)) # 178 out of 178 DArTag'ed accessions in pedigree
tagRelMat$stock_uniquename[!tagAccInPed]
#### DONE curating names

# This might not be important any more
# Fix indCallRate names
sameName <- names(indCallRate) %in% tagRelMat$stock_uniquename
print(sum(sameName))
names(indCallRate)[!sameName]
# "SA18-CB-S10-FG6" "SL20-JL-S3-FG1" "SL18-OI-S15-FG3" "SL18-UCONN-S142-MG2"
# "SL18-JS-S6-FG2" "SL18-NC-S5-FG3"
# I'm not sure why these ended up not being loaded.
# These are the first six samples of DArTag:
# SA18-CB-10-FG6,SL20-JL-3-FG-1,SL18-OI-15-FG3,SL18-UCONN-S142-MG2,SL18-JS-6-FG2,SL18-NC-5-FG3
# SL18-UCONN-S142-MG2 has a very bad call rate, but the others don't
# We should just go ahead and recalculate the tagRelMat using the markers
# Calculate GRM as (W %*% W^T) / sum(locusDosageVariance)
# To calculate p, the function expects dosage coding to be 0, 1, 2 for diploid
# and 0, 1 for haploid
### Calculate the relationship matrix for GPs
calcGenomicRelationshipMatrix <- function(locusMat, ploidy=2){
  if (!any(ploidy == 1:2)) stop("Ploidy must be 1 or 2")
  freq <- colMeans(locusMat, na.rm=T) / ploidy
  locusMat <- scale(locusMat, center=T, scale=F)
  return(tcrossprod(locusMat) / sum(ploidy*freq*(1-freq)))
}

# Deal with duplicated names.  NOTE: the duplicated names will not be
# in the pedigree. QualCtrlDArTag should have gotten rid of these
duplicated_names <- duplicated(colnames(dartagMrk), fromLast = FALSE)
if (any(duplicated_names)){
dartagMrk_DupCName <- dartagMrk %>%
  rename("SL20-MB-S2-FG3_2" = 178, "SL19-UCONN-S90-FG1_2" = 187)
colnames(dartagMrk_DupCName)
} else{
  dartagMrk_DupCName <- dartagMrk
}
# Recode DArTag markers
# Now they are 0 ref allele, 1 alt allele, 2 if both alleles.
# I am going to consider both alleles to be heterozygote, so recode to 0.5
# look into across() function
tst <- dartagMrk_DupCName %>%
  mutate_at(c(17:ncol(dartagMrk_DupCName)), funs(replace(., . > 1, 0.5))) # funs deprecated

source(here::here("analysis", "imputeMarkerMatrix.R"))

dartagMrkImputed <- impute.glmnet(t(tst[,17:ncol(tst)]))
hist(dartagMrkImputed[is.na(t(as.matrix(tst[,17:ncol(tst)])))])
colnames(dartagMrkImputed) <- tst$MarkerName

tagRelMat2 <- calcGenomicRelationshipMatrix(t(dartagMrkImputed), ploidy=1)

mrkToImpute <- t(tst[,17:ncol(tst)])
colnames(mrkToImpute) <- tst$MarkerName
Amat_rrBLUP <- rrBLUP::A.mat(as.matrix(mrkToImpute)*2 - 1, return.imputed=T,
                             impute.method="EM")
dartagMrkImpEM <- (Amat_rrBLUP$imputed + 1) / 2
dartagMrkImpEM[dartagMrkImpEM < 0] <- 0
dartagMrkImpEM[dartagMrkImpEM > 1] <- 1
toCompare <- tst %>% filter(MarkerName %in% colnames(dartagMrkImpEM))
hist(dartagMrkImpEM[is.na(t(as.matrix(toCompare[,17:ncol(toCompare)])))])
toCompare_glmnet <- dartagMrkImputed[, colnames(dartagMrkImputed) %in% colnames(dartagMrkImpEM)]

plot(toCompare_glmnet[is.na(t(as.matrix(toCompare[,17:ncol(toCompare)])))],
     dartagMrkImpEM[is.na(t(as.matrix(toCompare[,17:ncol(toCompare)])))], pch=16, cex=0.2)

cor(toCompare_glmnet[is.na(t(as.matrix(toCompare[,17:ncol(toCompare)])))],
    dartagMrkImpEM[is.na(t(as.matrix(toCompare[,17:ncol(toCompare)])))]) # 0.80

# Final DArTag relationship matrix here
dartImputed <- (toCompare_glmnet + dartagMrkImpEM) / 2
dartagRelMat <- calcGenomicRelationshipMatrix(dartImputed, ploidy=1)

sameName <- tagRelMat$stock_uniquename %in% names(indCallRate)
print(sum(sameName))
tagRelMat$stock_uniquename[!sameName]

# Final WGS relationship matrix here
wgsRelMat <- read_tsv(file=here::here("data", "WGSRelationshipMatrix.tsv"), na=c("", "NA", "N/A"))
print(sum(wgsRelMat$stock_uniquename %in% pedData$Accession))
# All the wgs accessions in pedigree

dim(fndRelMat); head(colnames(fndRelMat)); class(fndRelMat)
dim(dartagRelMat);  head(colnames(dartagRelMat)); class(dartagRelMat)
dim(wgsRelMat); head(colnames(wgsRelMat)); class(wgsRelMat)

# Make the pedigree relationship matrix
# Decide what to keep from the full pedigree data
# Keep:
# SPs for which we have phenotypes
# GPs for which we have marker data
# SPs that are parents of GPs for which we have marker data
# GPs that were used in crossing
# SPs that are parents of GPs that were used in crossing
# GPs that are progeny of SPs with phenotypes
# Marker datasets on GPs
gpsWithMarkers <- union(rownames(dartagRelMat), wgsRelMat$stock_uniquename)
spsWithMarkers <- fndRelMat$stock_uniquename

# Phenodata from GOM with decent photo scores
spsWithPhenotypes <- phenoData123G %>%
  filter(`Photo score 0-3|CO_360:0000300` > 0) %>% pull(germplasmName)

# Get SPs that are parents of GPs with marker data
pedGPsWithMrk <- pedData %>% filter(Accession %in% gpsWithMarkers)
parentsGPsWithMrk <- union(pedGPsWithMrk$Female_Parent, pedGPsWithMrk$Male_Parent) %>% setdiff(NA)

# Get GPs that are parents of SPs with phenotypes
pedSPsWithPhen <- pedData %>% filter(Accession %in% spsWithPhenotypes)
gpsParSPsWithPhen <- union(pedSPsWithPhen$Female_Parent, pedSPsWithPhen$Male_Parent) %>% setdiff(NA)

# Get SPs that are parents of those GPs
pedGPsParSPsWithPhen <- pedData %>% filter(Accession %in% gpsParSPsWithPhen)
grandParSPsWithPhen <- union(pedGPsParSPsWithPhen$Female_Parent, pedGPsParSPsWithPhen$Male_Parent) %>% setdiff(NA)

# Get GPs that are progeny of SPs with phenotypes
gpsProgOfSPsWithPhen <- pedData %>%
  filter(Female_Parent %in% spsWithPhenotypes |
           Male_Parent %in% spsWithPhenotypes) %>%
  pull(Accession) %>% setdiff(NA)

allGPsToKeep <- union(gpsWithMarkers, gpsParSPsWithPhen) %>%
  union(gpsProgOfSPsWithPhen)
allSPsToKeep <- union(spsWithMarkers, parentsGPsWithMrk) %>%
  union(grandParSPsWithPhen) %>%
  union(spsWithPhenotypes)

pedKeep <- pedData %>% dplyr::filter(Accession %in% allGPsToKeep |
                                       Accession %in% allSPsToKeep)

# To be ahead of its progeny, SL-ME-MG1 has to go higher up
meMG1 <- pedKeep %>% dplyr::filter(Accession == "SL-ME-MG1")
pedKeep <- pedKeep %>% dplyr::filter(Accession != "SL-ME-MG1")
pedKeep <- pedKeep %>% add_row(meMG1, .after=1)

# Try to assign ploidy level
pedKeep <- pedKeep %>%
  dplyr::mutate(ploidy=if_else(Accession %in% allGPsToKeep, "GP", "SP"))
# Do some verifying
# Look at GP rows that do not have "FG" or "MG"
fgOrMG <- c(grep("FG", pedKeep$Accession), grep("MG", pedKeep$Accession))
noFGMG <- dplyr::setdiff(1:nrow(pedKeep), fgOrMG)
gpNoFGMG <- dplyr::intersect(which(pedKeep$ploidy == "GP"), noFGMG)
# There are SPs with names containing "-Mixed" labeled as GPs
pedKeep$ploidy[grep("-Mixed", pedKeep$Accession, fixed=T)] <-  "SP"

# View(pedKeep %>% filter(ploidy=="GP"))
# What are these?: SL18-CT1-SFG-3, SL18-CT1-SMG-2, SL18-CT1-SMG-3
# Missing parents: SL18-FW-S12-MG11, SL18-ME-S1-MG2, SL18-OI-S15-FEMALE
# SL-ME-MG1
# View(pedKeep %>% filter(ploidy=="SP"))
# Missing parents: SL19-UCONN-107, SL19-UCONN-158, SL19-UCONN-176
# What are these?: SL19-UCONN-C1, SL19-UCONN-C2, SL21-UCONN-C1, SL21-UCONN-C2

# Let's try convertNamesToRows
source(here::here("code", "calcRelationshipMatrices.R"))
pedCNR <- pedKeep
pedCNR <- convertNamesToRows(as.matrix(pedCNR))
# Check whether there are parents that come after their offspring
rowNum <- 1:nrow(pedCNR)
sum(pedCNR[, 2] > rowNum)
sum(pedCNR[, 3] > rowNum)

source(here::here("code", "GS_ms", "calcCCmatrixBiphasic.R"))
pedCNR[pedKeep$ploidy == "GP", 3] <- NA
pedRelMat <- calcCCmatrixHaploid(pedCNR)

if (FALSE){ # I don't think this is strictly needed. To move forward skip for now
  # I want to order the pedigree a little
  # SPs in 2019 to 2022 order
  # Add time for GPs that were offspring of SPs
  # Substract time for GPs that were parents of SPs
  # Put all the SP founders at the beginning
  # iy is individual year
  iy <- pedKeep$Accession %>% strsplit("-", fixed=T) %>%
    sapply(function(sv) return(as.numeric(substr(sv[1], 3, 4))))
  iy[is.na(iy)] <- 0
  pedKeep <- pedKeep %>% dplyr::mutate(indRelTime=iy*10)

  # Get GPs that are progeny of SPs with phenotypes
  gpsProgOfSPsWithPhen <- pedData %>%
    filter(Female_Parent %in% spsWithPhenotypes |
             Male_Parent %in% spsWithPhenotypes) %>%
    pull(Accession) %>% setdiff(NA)

  # Get SPs that are parents of GPs with marker data
  pedGPsWithMrk <- pedData %>% filter(Accession %in% gpsWithMarkers)
  parentsGPsWithMrk <- union(pedGPsWithMrk$Female_Parent, pedGPsWithMrk$Male_Parent) %>% setdiff(NA)

  # Get GPs that are parents of SPs with phenotypes
  pedSPsWithPhen <- pedData %>% filter(Accession %in% spsWithPhenotypes)
  gpsParSPsWithPhen <- union(pedSPsWithPhen$Female_Parent, pedSPsWithPhen$Male_Parent) %>% setdiff(NA)

  # Get SPs that are parents of those GPs
  pedGPsParSPsWithPhen <- pedData %>% filter(Accession %in% gpsParSPsWithPhen)
  grandParSPsWithPhen <- union(pedGPsParSPsWithPhen$Female_Parent, pedGPsParSPsWithPhen$Male_Parent) %>% setdiff(NA)
}

###############################################################################
# Work to combine the covariance matrices here
wgsRelMatm <- as.matrix(wgsRelMat[,-1])
rownames(wgsRelMatm) <- colnames(wgsRelMatm)
fndRelMatm <- as.matrix(fndRelMat[,-1])
rownames(fndRelMatm) <- colnames(fndRelMatm)
pedRelMatm <- pedRelMat$mccMat
rownames(pedRelMatm) <- colnames(pedRelMatm) <- pedKeep$Accession
# Some checks to make sure this is all good
all(rownames(dartagRelMat)==colnames(dartagRelMat))
all(rownames(pedRelMatm)==colnames(pedRelMatm))
all(rownames(dartagRelMat) %in% rownames(pedRelMatm))
rownames(dartagRelMat)[!(rownames(dartagRelMat) %in% rownames(pedRelMatm))]
all(rownames(fndRelMatm) %in% rownames(pedRelMatm))
all(rownames(wgsRelMatm) %in% rownames(pedRelMatm))

# The pedigree-based matrix is a coefficient of coancestry matrix
# The fndRelMat came from Sugarkelpbase and is on diploids. Because I want to
# combine coefficients of coancestry, I need to divide it by two
fndRelMatm <- fndRelMatm/2
# The dartagRelMat was calculated from DArTag markers, specifying ploidy=1. This
# one is all set
# The wgsRelMat came from Sugarkelpbase and is on haploids. Because this was
# calculated assuming the genotypes were diploid, I need to divide it by two
wgsRelMatm <- wgsRelMatm/2

# Function to return offdiagonal elements of a relationship matrix
offDiag <- function(relMat){
  return(c(relMat[upper.tri(relMat)]))
}

plotTwoRelMat <- function(relMat1, relMat2){
  if (nrow(relMat1) != ncol(relMat1)) stop("Matrices have to be square")
  if (nrow(relMat1) != nrow(relMat2)) stop("Matrices have to have same size")
  if (ncol(relMat1) != ncol(relMat2)) stop("Matrices have to have same size")
  plot(offDiag(relMat1), offDiag(relMat2), pch=16,
       main="relMat1 against relMat2, off-diagonal",
       xlab="relMat1 pairwise", ylab="relMat2 pairwise")
}

plotTwoRelMat(pedRelMatm[rownames(fndRelMatm), colnames(fndRelMatm)],
              fndRelMatm)
plotTwoRelMat(pedRelMatm[rownames(wgsRelMatm), colnames(wgsRelMatm)],
              wgsRelMatm)
plotTwoRelMat(pedRelMatm[rownames(dartagRelMat), colnames(dartagRelMat)],
              dartagRelMat)

hist(diag(fndRelMatm))
hist(diag(wgsRelMatm))
hist(diag(dartagRelMat))
hist(diag(pedRelMatm[rownames(fndRelMatm), colnames(fndRelMatm)]))
hist(diag(pedRelMatm[rownames(dartagRelMat), colnames(dartagRelMat)]))
unique(diag(pedRelMatm[pedKeep$ploidy=="SP", pedKeep$ploidy=="SP"]))
unique(diag(pedRelMatm[pedKeep$ploidy=="GP", pedKeep$ploidy=="GP"]))

# OK. None of this is terribly conclusive that we've done it all right, but it's
# what we've got
library(CovCombR)

covList<-NULL
covList[[1]] <- pedRelMatm + diag(1e-3, nrow(pedRelMatm))
covList[[2]] <- fndRelMatm + diag(1e-3, nrow(fndRelMatm))
covList[[3]] <- wgsRelMatm + diag(1e-3, nrow(wgsRelMatm))
covList[[4]] <- dartagRelMat + diag(1e-3, nrow(dartagRelMat))

weights <- c(1, 2, 2, 2)
# The CovComb part needs to be done on a server with many cores...
# outCovComb <- CovComb(Klist=covList, nu=3000, w=weights, Kinit=covList[[1]])
# saveRDS(outCovComb, here::here("output", "outCovComb.rds"))

outCovComb <- readRDS(here::here("output", "outCovComb.rds"))

# Looking at a few parameters to make sure they make sense
phenoData123G$studyYear %>% unique
phenoData123G$studyName %>% unique
phenoData123G$locationName %>% unique
phenoData123G$germplasmName %>% unique %>% length
(phenoData123G$germplasmName %>% unique) %in% (outCovComb %>% rownames) %>% sum

# Linear model:
# y = Xb + Zu
# The fixed effects will only by the trials ($studyName)
X <- model.matrix(as.formula("~ -1 + studyName"), data=phenoData123G)
# Make the incidence matrix connecting outCovComb with phenoData123G
allLevFact <- rownames(outCovComb)
Z <- model.matrix(~ -1 + factor(germplasmName, levels=allLevFact),
                  data=phenoData123G)

# Proof of concept
tst <- rrBLUP::mixed.solve(y=phenoData123G$`Wet yield kg/m|CO_360:0000304`,
                           Z=Z, K=outCovComb, X=X)
head(tst$u)

# Generalize
phenoData123G <- phenoData123G %>%
  dplyr::mutate(logDYld=log(`Dry yield kg/m|CO_360:0000310` + 1))
phenoData123G <-  phenoData123G %>%
  dplyr::mutate(logWYld=log(`Wet yield kg/m|CO_360:0000304` + 1))
notIsNAphenDat <- apply(phenoData123G, 2, function(v) sum(!is.na(v)))
traitsToAnalyze <- phenoData123G %>%
  dplyr::select((notIsNAphenDat > 300) %>% which) %>%
  dplyr::select(contains("CO_360")) %>% colnames
traitsToAnalyze <- c(traitsToAnalyze, c("logDYld", "logWYld"))

calcBLUPsOnTrait <- function(trait){
  print(trait)
  y <- pull(phenoData123G, trait)
  print(range(y, na.rm=T))
  print(paste("Number equal to Zero", sum(y == 0, na.rm=T)))
  rrBLUPoutTrait <- try(
    rrBLUP::mixed.solve(y=y, Z=Z, K=outCovComb, X=X),
    silent=TRUE)
  if(inherits(rrBLUPoutTrait, "try-error")){
    toRet <- numeric(nrow(outCovComb))
  } else{
    toRet <- rrBLUPoutTrait$u
  }
  return(toRet)
}

allBLUPs <- sapply(traitsToAnalyze, calcBLUPsOnTrait)

# apply(allBLUPs, 2, range)
tst <- cor(allBLUPs)
View(tst)

allBLUPst <- bind_cols(pedKeep,
                       as_tibble(allBLUPs))
readr::write_csv(allBLUPst,
                 here::here("output", "allBLUPs.csv"),
                 quote="none", col_names=T)

###############################################################################
# There are some accessions where the SP S number has leading zero, which
# apparently I don't like.  They also exist without the leading zero.
curatePed <- FALSE
if (curatePed){
  leadZ <- pedData$Accession %>% sapply(function(s) gregexpr("S0", s, fixed=T)[[1]] != -1)
  lz <- pedData %>% filter(leadZ)
  leadZmp <- pedData$Male_Parent %>% sapply(function(s) gregexpr("S0", s, fixed=T)[[1]] != -1)
  lzmpt <- pedData %>% filter(leadZmp)
  leadZfp <- pedData$Female_Parent %>% sapply(function(s) gregexpr("S0", s, fixed=T)[[1]] != -1)
  lzfpt <- pedData %>% filter(leadZfp)

  leadZfp <- pedData$Female_Parent == "SL18-CC-S1"
  for (i in 2:9){
    leadZfp <- leadZfp | pedData$Female_Parent == paste0("SL18-CC-S", i)
  }
  lzfpt <- pedData %>% filter(leadZfp)

  leadZmp <- pedData$Male_Parent == "SL18-CC-S1"
  for (i in 2:9){
    leadZmp <- leadZmp | pedData$Male_Parent == paste0("SL18-CC-S", i)
  }
  lzmpt <- pedData %>% filter(leadZmp)

  ccS <- pedData$Accession %>% sapply(function(s) gregexpr("SL18-CC-S", s, fixed=T)[[1]] != -1)
  ccSt <- pedData %>% filter(ccS)
  nDash <- ccSt$Accession %>% sapply(function(s) gregexpr("-", s, fixed=T)[[1]] %>% length)
  ccStSP <- ccSt %>% filter(nDash == 2)
}
