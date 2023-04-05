#################################################
## This function uses glmnet to impute NA values
#################################################

library(glmnet)

# NOTES
# You enter the function with a matrix, individuals in rows, markers in columns
# I have only tested it with markers coded as -1, 0, 1 for AA, AB, and BB, but
# other codings would probably work too.
# Returns a matrix of the same dimensions, but with no missing data. Imputed values
# are real numbers (not integers). This may be problematic for downstream mapping software.
# by default glmnet will look at ~100 different lambda penalty coefficients.
# It approximately doubles the speed to look at only 10 values. That probably
# lowers the accuracy by a couple percent, but not much more.
# Another thing that ~ doubles the speed is to do  5-fold rather than 10-fold cv
# The thing that makes the most difference is not putting ALL the other
# markers in as predictors, but only the top xx of them. I am using 60 now.
impute.glmnet <- function(matNA){
  cvLambda <- exp(-(2:11))
  # Start with mean impute
  matNoNA <- apply(matNA, 2, function(vec){vec[is.na(vec)] <- mean(vec, na.rm=TRUE); return(vec)})
  # I am using 60 markers for prediction.  We could experiment with this parameter.
  # Might crash if you have fewer than nPred markers in the matrix.
  nPred <- min(60, round(ncol(matNA) * 0.5))
  # Function to fill in one column with glmnet impute
  imputeMrk <- function(k){
    varRange <- range(matNA[,k], na.rm=TRUE) # Use to prevent imputations from going outside the original range
    isNA <- is.na(matNA[,k])
    # If the marker is monomorphic, impute with the sole value
    if (sd(matNA[,k], na.rm=TRUE) == 0) matNoNA[isNA,k] <<- matNA[which(!isNA)[1],k] else{
      corMrk <- abs(cov(matNA[,k], matNA, use="pairwise.complete.obs"))
      # Retain markers that correlate highly with marker to be imputed
      predMrk <- setdiff(order(corMrk,decreasing=TRUE)[1:nPred], k)
      cvModels <- cv.glmnet(x=matNoNA[!isNA,predMrk], y=matNA[!isNA,k], nfolds=5, lambda=cvLambda)
      # The double assignment arrow puts values into matNoNA defined above
      pred <- predict(cvModels, s="lambda.min", newx=matNoNA[isNA,predMrk, drop=FALSE])
      pred[pred < varRange[1]] <- varRange[1]
      pred[pred > varRange[2]] <- varRange[2]
      matNoNA[isNA,k] <<- pred
    }
    return(k)
  }
  # Go in order from least to most missing (probably not really needed)
  sumIsNA <- apply(matNA, 2, function(vec) sum(is.na(vec)))
  imputeOrder <- order(sumIsNA)
  imputeOrder <- imputeOrder[sumIsNA[imputeOrder] > 0] # Don't impute if none missing
  dummy <- try(sapply(imputeOrder, imputeMrk), silent=TRUE)
  return(matNoNA)
}

# Create nImpute sets of matrix cells to be masked
# Each set should be about the same size and no set should have the same cell twice
# This is pretty tricky.  I wonder if there is a better way to do it.
setupMaskGroupList <- function(notNA, nImpute, nImputeNotNA){
  nNotNA <- length(notNA)
  totImpute <- nImputeNotNA*nNotNA
  groupSizes <- rep(totImpute %/% nImpute, nImpute)
  addOne <- sample(nImpute, totImpute %% nImpute)
  groupSizes[addOne] <- groupSizes[addOne] + 1
  currIdx <- 1:nNotNA
  maskGroupList <- list()
  for (group in 1:nImpute){
    thisGroup <- NULL
    nCurrIdx <- length(currIdx)
    if (nCurrIdx == groupSizes[group]){
      thisGroup <- currIdx
      currIdx <- 1:nNotNA
    } else{
      if (nCurrIdx < groupSizes[group]){
        thisGroup <- currIdx
        currIdx <- 1:nNotNA
        sampleFrom <- currIdx[-thisGroup]
        sampleSize <- groupSizes[group] - nCurrIdx
      } else{
        sampleFrom <- nCurrIdx
        sampleSize <- groupSizes[group]
      }
      thisSample <- sample(sampleFrom, sampleSize)
      thisGroup <- c(thisGroup, currIdx[thisSample])
      currIdx <- currIdx[-thisSample]
    }
    maskGroupList <- c(maskGroupList, list(thisGroup))
  }
  return(maskGroupList)
}
# This function uses impute.glmnet to reimpute data that are _not_ missing and potentially identify
# observed values that are not supported by the remaining values in the matrix
# Reasonable values of nImpute and nImputeNotNA are probably something like 100 and 5
# Those values mean that each observed value (notNA) will be reimputed 5 times over the course of imputing
# the complete dataset 100 times.  Imputing the complete dataset that many times will obviously take a while.
# The fraction of observed values that will be masked in any given imputation is nImputeNotNA / nImpute.
imputationAndErrorFlagging <- function(matNA, nImpute, nImputeNotNA){
  isNA <- which(is.na(matNA))
  notNA <- 1:length(matNA)
  if (length(isNA) > 0) notNA <- notNA[-isNA]
  maskGroupList <- setupMaskGroupList(notNA, nImpute, nImputeNotNA)

  imputeWithMasking <- function(maskIdx){
    matNA[maskIdx] <- NA # This does not change matNA outside of this function
    matNoNA <- impute.glmnet(matNA)
    return(list(imputedNA=matNoNA[isNA], imputedNotNA=matNoNA[maskIdx]))
  }

  nNotNA <- length(notNA)
  imputedNA <- matrix(nrow=length(isNA), ncol=nImpute)
  imputedNotNA <- matrix(nrow=nNotNA, ncol=nImputeNotNA)
  whichRep <- rep(1, nNotNA)
  for (group in 1:nImpute){
    cat("Doing group", group, "\n")
    maskIdx <- maskGroupList[[group]]
    imputeOut <- imputeWithMasking(notNA[maskIdx])
    imputedNA[,group] <- imputeOut$imputedNA
    imputedNotNA[cbind(maskIdx, whichRep[maskIdx])] <- imputeOut$imputedNotNA
    whichRep[maskIdx] <- whichRep[maskIdx] + 1
  }

  matNA[isNA] <- rowMeans(imputedNA)
  diffObsImp <- matNA[notNA] - rowMeans(imputedNotNA)
  return(list(isNA=isNA, notNA=notNA, matNoNA=matNA, diffObsImp=diffObsImp, iNA=imputedNA, iNNA=imputedNotNA))
}

# Multicore (parallel) version of the imputation method above
# NOTE: this cannot be run from the R or RStudio GUI.  Has to be run from command line
# NOTE: I have not tested this multicore version yet (23 Apr 2014).  Should work though.
imputationAndErrorFlagging.mc <- function(matNA, nImpute, nImputeNotNA, nCores=4){
  library("multicore")

  isNA <- which(is.na(matNA))
  notNA <- 1:length(matNA)
  if (length(isNA) > 0) notNA <- notNA[-isNA]
  maskGroupList <- setupMaskGroupList(notNA, nImpute, nImputeNotNA)

  imputeWithMasking <- function(maskIdx){
    maskIdx <- notNA[maskIdx]
    matNA[maskIdx] <- NA # This does not change matNA outside of this function
    matNoNA <- impute.glmnet(matNA)
    return(list(imputedNA=matNoNA[isNA], imputedNotNA=matNoNA[maskIdx]))
  }
  # Do the imputations
  imputeOut <- mclapply(maskGroupList, imputeWithMasking, mc.cores=nCores)

  # Collect the output into matrices
  nNotNA <- length(notNA)
  imputedNA <- matrix(nrow=length(isNA), ncol=nImpute)
  imputedNotNA <- matrix(nrow=nNotNA, ncol=nImputeNotNA)
  whichRep <- rep(1, nNotNA)
  for (group in 1:nImpute){
    maskIdx <- maskGroupList[[group]]
    imputedNA[,group] <- imputeOut[[group]]$imputedNA
    imputedNotNA[cbind(maskIdx, whichRep[maskIdx])] <- imputeOut[[group]]$imputedNotNA
    whichRep[maskIdx] <- whichRep[maskIdx] + 1
  }

  matNA[isNA] <- rowMeans(imputedNA)
  diffObsImp <- matNA[notNA] - rowMeans(imputedNotNA)
  return(list(isNA=isNA, notNA=notNA, matNoNA=matNA, diffObsImp=diffObsImp, iNA=imputedNA, iNNA=imputedNotNA, imputeOut=imputeOut))
}

#############################################
# Start here: scripts using the functions
#############################################
biparental <- FALSE
if (biparental){
  # Run simple tests on a small part of the biparental snp scores
  load("cassavaBiparFam.RData") # You should use the "stringent" and "lax" calls from the vcf that you set up in HW1
  testMat <- t(snpsScores[1:100, 1:50])
} else{
  # Use this to try out the system on a panel population rather than a biparental family
  load("./IITA cassava data/IITA - SNP Data - Cycle 0 - ggayt.RData")
  testMat <- ggaytNA[1:50, 1:100]
}
# Find out how long it takes for this small matrix (4 minutes on my computer)
s <- Sys.time()
testOut <- imputationAndErrorFlagging(testMat, 40, 2)
print(Sys.time() - s)
# Identify on this small matrix which cells have high discordance
# That is, where the imputation is far from the observed
hasHighDiscord1 <- testOut$notNA[which(abs(testOut$diffObsImp) > 1)]
cat("Cells in the matrix with high discordance", hasHighDiscord1, "\n")
cat("Num not NA with high discordance:", length(hasHighDiscord1), "\n")
# Handier to know the row and column in the matrix rather then the vector position
rcHiDiscord <- (hasHighDiscord1 - 1) %% 50 + 1
rcHiDiscord <- cbind(rcHiDiscord, (hasHighDiscord1 - 1) %/% 50 + 1)
colnames(rcHiDiscord) <- c("row", "col")

# If you take a marker data set and spike in errors, will the error flagging routine find them?
nErr <- 20 # Number of errors to spike in
testResult <- NULL
for (repetition in 1:3){
  spikeMat <- testMat
  # Spike in nErr big errors by flipping 0 to 2 and vice versa
  is0or2 <- which(testMat == 0 | testMat == 2)
  spiked <- sample(is0or2, nErr)
  spikeMat[spiked] <- 2 - spikeMat[spiked]
  # Run the function to flag potential observation errors
  testOut <- imputationAndErrorFlagging(spikeMat, 40, 2)
  # Find out how many errors were found based on the difference between the imputed and observed values
  stdDevDiscord <- sd(testOut$diffObsImp, na.rm=TRUE)
  # High imputation error is arbitrarily diffObsImp > 1
  hasHighDiscord <- testOut$notNA[abs(testOut$diffObsImp) > 1]
  cat("Num not NA with high discordance:", length(hasHighDiscord), "\n")
  cat("Num errors found out of", nErr, ":", sum(spiked %in% hasHighDiscord), "\n")
  testResult <- rbind(testResult, c(stdDevDiscord, length(hasHighDiscord), sum(spiked %in% hasHighDiscord)))
}

# All these flips were found every time (I think)
colnames(testResult) <- c("stdDevDiscord", "nHighDiscordance", "nTrueErrorAmongHighDiscordance")
print(colMeans(testResult))
