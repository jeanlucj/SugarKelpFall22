library(tidyverse)
here::i_am("analysis/QualCtrlDArTag.R")

colTypes <- paste0(paste0(rep("c", 4), collapse=""),
                   paste0(rep("d", 12), collapse=""),
                   paste0(rep("i", 188), collapse=""))
# As best I can tell, 0 = Hom Ref, 1 = Hom Alt, 2 = Het
dartagMrk <- read_csv(file=here::here("data", "Report_DSacc22-7080_SNP.csv"),
                skip = 6, na = "-", col_types=colTypes)

# Some curating of column names
cnames <- colnames(dartagMrk)
cnames <- gsub("Old", "OLD", cnames, fixed=T)
cnames <- gsub("G-", "G", cnames, fixed=T)
cnames <- gsub("UCONN-S", "UCONN-", cnames, fixed=T)
cnames <- gsub("UCONN-", "UCONN-S", cnames, fixed=T)
cnames[cnames == "SA-CB-7-MG3"] <- "SA18-CB-7-MG3"
missingS <- substr(cnames, 8, 8) == "-"
addS <- paste0(substr(cnames[missingS], 1, 8), "S", substr(cnames[missingS], 9, 100))
cnames[missingS] <- addS
cnames <- gsub("SFG", "FG", cnames)
colnames(dartagMrk) <- cnames
# Remove some data that corresponds to GPs with ambiguity
removeFromTagData <- which(cnames %in% c("SL19-UCONN-S90-FG1-2", "SL20-MB-S11-FGOLD1", "SL20-MB-S2-FG3-2", "SL20-MB-S6-MGOLD12"))
dartagMrk <- dartagMrk[,-removeFromTagData]
# NOTE: We might want to change what we do here...
# Remove duplicated clones
dartagMrk <- dartagMrk[,!duplicated(colnames(dartagMrk))]

mrkCallRate <- apply(dartagMrk[, 17:ncol(dartagMrk)], 1, function(v) sum(!is.na(v))/188)
jpeg("./output/HistMrkCallRate.jpeg")
hist(mrkCallRate)
dev.off()

indCallRate <- apply(dartagMrk[, 17:ncol(dartagMrk)], 2, function(v) sum(!is.na(v))/nrow(dartagMrk))
jpeg("./output/IndMrkCallRate.jpeg")
hist(indCallRate)
dev.off()

badIndCall <- indCallRate < 0.1
bic <- tibble(gp_name=names(indCallRate)[badIndCall], call_rate=indCallRate[badIndCall]) %>%  write_csv(file=here::here("output", "IndLowCallRate.csv"))

mrkFreqHet <- apply(dartagMrk[, 17:ncol(dartagMrk)], 1, function(v) sum(v == 2, na.rm=T) / sum(!is.na(v)))
hist(mrkFreqHet)
plot(dartagMrk$AvgPIC, mrkFreqHet)

indFreqHet <- apply(dartagMrk[, 17:ncol(dartagMrk)], 2, function(v) sum(v == 2, na.rm=T) / sum(!is.na(v)))
jpeg("./output/IndFreqHet.jpeg")
hist(indFreqHet)
dev.off()

badIndHet <- indFreqHet > 0.03
bih <- tibble(gp_name=names(indFreqHet)[badIndHet], het_rate=indFreqHet[badIndHet]) %>%
  write_csv(file="./output/IndHighHetRate.csv")

sum(indFreqHet > 0.03)

jpeg("./output/HetRateAgainstCallRate.jpeg")
plot(indCallRate, indFreqHet)
dev.off()
