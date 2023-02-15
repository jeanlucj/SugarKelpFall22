setwd("~/Documents/GitRepo/SugarKelpFall22")
here::i_am("analysis/pullInData.R")
library(tidyverse)

source(here::here("analysis", "QualCtrlDArTag.R"))
# Phenotypes
# This file came from downloading all phenotypes from trials where the
# name includes "FARM"
phenoData <- read_csv(file=here::here("data", "AllFarmPhenotypes.csv"),
                      skip=3, guess_max=10000, na=c("", "NA", "N/A"))
phenoData <- phenoData %>% filter(observationLevel == "plot")
sumIsNA <- apply(phenoData, 2, function(v) sum(is.na(v)))
phenoData <- phenoData %>% select(which(sumIsNA<nrow(phenoData)))
phenoData <- phenoData %>% select(-contains("DbId"))
phenoData <- phenoData %>% select(-contains("program"))
justTraits <- phenoData %>% select(contains("CO_360")) %>%
  select(-contains(c("253", "320")))
hasData <- apply(justTraits, 1, function(v) sum(!is.na(v))) > 2
phenoData <- phenoData %>% filter(hasData)

# Pedigree data
pedData <- read_tsv(file=here::here("data", "AllAccessionsPedigree.txt"),
                    na=c("", "NA", "N/A"))
# Manual curating
pedData <- pedData %>% filter(Accession != "SA18-CB-FG5")
nDash <- pedData$Accession %>% sapply(function(s) gregexpr("-", s, fixed=T)[[1]] %>% length)
tstSP <- pedData %>% filter(nDash == 2)
tstGP <- pedData %>% filter(nDash > 2)

# NOTE: I have to do some curation using the markers on the DArTag to eliminate
# the ones that look like they are mixtures

# Marker relationship matrices
fndRelMat <- read_tsv(file=here::here("data", "DArTSeqRelationshipMatrix029910.tsv"), na=c("", "NA", "N/A"))
fndAccInPed <- fndRelMat$stock_uniquename %in% pedData$Accession
print(sum(fndAccInPed)) # All 179 founders in pedigree

tagRelMat <- read_tsv(file=here::here("data", "DArTagRelationshipMatrix029910.tsv"), na=c("", "NA", "N/A"))
tagAccInPed <- tagRelMat$stock_uniquename %in% pedData$Accession
print(sum(tagAccInPed)) # 172 out of 182 DArTag'ed accessions in pedigree
tagRelMat$stock_uniquename[!tagAccInPed]
# [1] "SL18-SF-S13-SFG4"     "SL18-SF-S13-SFG5"     "SL18-SF-S13-SFG6"
# [4] "SL19-UCONN-S156-FG3"  "SL19-UCONN-S158-FG3"  "SL19-UCONN-S23-FG4"
# [7] "SL19-UCONN-S90-FG1-2" "SL20-MB-S11-FGOLD1"   "SL20-MB-S2-FG3-2"
#[10] "SL20-MB-S6-MGOLD12"
# Remove the SFG names
tagRelMat$stock_uniquename <- gsub("SFG", "FG", tagRelMat$stock_uniquename)
# These "SL19-UCONN-S156-FG3"  "SL19-UCONN-S158-FG3"  "SL19-UCONN-S23-FG4"
# Have exact matches in the database: not sure why they are not in the pedigree

# Fix indCallRate names
sameName <- names(indCallRate) %in% tagRelMat$stock_uniquename
print(sum(sameName))
names(indCallRate)[!sameName]
# "SA18-CB-S10-FG6"     "SL20-JL-S3-FG1"      "SL18-OI-S15-FG3"     "SL18-UCONN-S142-MG2"
# "SL18-JS-S6-FG2"      "SL18-NC-S5-FG3"
# I'm not sure why these ended up not being loaded.
# SL18-UCONN-S142-MG2 has a very bad call rate, but the others don't

sameName <- tagRelMat$stock_uniquename %in% names(indCallRate)
print(sum(sameName))
tagRelMat$stock_uniquename[!sameName]

wgsRelMat <- read_tsv(file=here::here("data", "WGSRelationshipMatrix.tsv"), na=c("", "NA", "N/A"))
wgsRelMat <- wgsRelMat$stock_uniquename %in% pedData$Accession
# All the wgs accessions in pedigree

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
}
