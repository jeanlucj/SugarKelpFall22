#' Calculate a coefficient of coancestry matrix
#'
#' \code{calcCCmatrixBiphasic} returns an additive relationship matrix from a
#'  pedigree specified in three columns. The first column has to be the row
#'  number and sire and dam columns refer directly to rows of the pedigree.
#'
#' \code{calcCCmatrixBiphasic} has functionality useful for species with a
#'  biphasic lifecycle. A haploid progeny can be specified by giving the
#'  row number of the diploid parent in the first column and setting the
#'  second column to NA. A diploid progeny should be specified by giving
#'  the row numbers of the two haploid parents.
#'
#'  The following rules are followed.
#'  1. Diploid from two diploids. The usual as given in quantitative genetics.
#'  2. Diploid from two haploids. The row is the mean of the haploid rows.
#'  The diagonal element is 0.5 * (1 + coancestry of the two haploids).
#'  3. Haploid from diploid. The row is the same as the row of the diploid.
#'  The diagonal element is 1.
#'  4. Haploid from haploid. The row is the same as the row of the haploid.
#'  The diagonal element is 1. (Equivalent to cloning the haploid).
#'
#' @param pedColumns A data.frame with three columns. The first column
#'  has to be the row number and sire and dam columns refer directly to rows
#'  of the pedigree. Parents of founders need to be set to 0 (ZERO). The row
#'  of a child has to be after (i.e. a higher row number) that of its parents.
#'  Indicate a haploid progeny by setting the second column to NA. If both
#'  columns have numbers in them, the progeny is assumed to be diploid.
#'  If an individual has one known and one unknown parent, set the unknown
#'  parent to 0.
#'
#' @return A matrix, \code{ccMat}, the coefficient of coancestry matrix
#'
calcCCmatrixBiphasic <- function(pedColumns){
    calcCCmatRow <- function(pedRec){ # Function to process one row of pedColumns
        prog <- pedRec[1]
        sire <- pedRec[2]
        dam <- pedRec[3]
        progM1 <- prog - 1
        if (sire != 0){
            sireRow <- ccMat[sire, 1:progM1] # Non-founder
        } else{
            sireRow <- integer(progM1) # Founder
        }
        if (is.na(dam)){ # The progeny is a haploid (GP)
            ccMat[prog, 1:progM1] <<- sireRow
            ccMat[1:progM1, prog] <<- sireRow
            ccSelf <- 1
        } else{ # The progeny is a diploid (SP)
            if (dam != 0){
                damRow <- ccMat[dam, 1:progM1] # Non-founder
            } else{
                damRow <- integer(progM1) # Founder
            }
            ccMat[prog, 1:progM1] <<- (sireRow + damRow) / 2
            ccMat[1:progM1, prog] <<- (sireRow + damRow) / 2
            if (sire > 0 & dam > 0){
                ccSelf <- (1 + ccMat[sire, dam]) / 2
            } else{
                ccSelf <- 0.5
            }
        }
        ccMat[prog, prog] <<- ccSelf
    }#END calcCCmatRow
    
    # calculate Coef Coan matrix here
    nInd <- nrow(pedColumns)
    ccMat <- matrix(0, nInd, nInd)
    ccMat[1, 1] <- ifelse(is.na(pedColumns[1, 3]), 1, 0.5)
    tmp <- apply(pedColumns[2:nInd,], 1, calcCCmatRow)
    rownames(ccMat) <- colnames(ccMat) <- pedColumns[,1]
    return(ccMat)
}

#' @param gMat A genomic relationship matrix.
#'
#' @param aMat A pedigree relationship matrix.
#'
#' @param aMatFounders The individuals to match in diagonal with gMat.
#'
#' @param relWgt Relative weights of the two matrices.
#'
#' @return A matrix, \code{ccMat}, the coefficient of coancestry matrix
#'
calcHmatrix <- function(gMat, aMat, aMatFounders, wgtAmat=0.05){
    # Scale the gMat so its mean diagonal is equal to the mean diagonal of aMatFounders
    amfMean <- mean(diag(aMat[aMatFounders, aMatFounders]))
    gMatMean <- mean(diag(gMat))
    gMat <- gMat * amfMean / gMatMean
    # If there are individuals in gMat but not aMat, they will be made founders
    gNoA <- setdiff(rownames(gMat), rownames(aMat))
    if (length(gNoA) > 0){
        nGnoA <- length(gNoA)
        nAmat <- nrow(aMat)
        aMat <- cbind(rbind(diag(amfMean, nGnoA), matrix(0, nAmat, nGnoA)), rbind(matrix(0, nGnoA, nAmat), aMat))
        rownames(aMat)[1:nGnoA] <- colnames(aMat)[1:nGnoA] <- gNoA
    }
    # Make sure things are set up right to create the partition
    idxGinA <- sapply(rownames(gMat), function(n) which(rownames(aMat) == n))
    newOrder <- c(idxGinA, setdiff(1:nrow(aMat), idxGinA))
    aMat <- aMat[newOrder, newOrder]
    # Do the math
    idx11 <- 1:nrow(gMat)
    a11 <- aMat[idx11, idx11]
    a11inva12 <- solve(a11) %*% aMat[1:nrow(gMat), -(1:nrow(gMat))]
    gw <- (1 - wgtAmat)*gMat + wgtAmat*a11
    hMat <- cbind(rbind(gw, crossprod(a11inva12, gw)),
    rbind(gw, crossprod(a11inva12, gw - a11)) %*% a11inva12)
    hMat[-idx11, -idx11] <- hMat[-idx11, -idx11] + aMat[-idx11, -idx11]
    rownames(hMat) <- colnames(hMat) <- rownames(aMat)
    return(hMat)
}

#' Calculate a haploid coefficient of coancestry matrix
#'
#' \code{calcCCmatrixHaploid} returns an additive relationship matrix from a
#'  pedigree specified in three columns. The first column has to be the row
#'  number and sire and dam columns refer directly to rows of the pedigree.
#'  The three-column pedigree can have haploids and diploids in it. It will
#'  be expanded to a two-column pedigree that just has haploids. In that
#'  pedigree, each diploid will be represented by two rows. The function will
#'  return this haploid relationship matrix.
#'  A haploid progeny can be specified by giving the row number of the diploid
#'  parent in the first column and setting the second column to NA. A diploid
#'  progeny should be specified by giving the row numbers of the two haploid
#'  parents.
#'
#'  The following rules are followed.
#'  1. When a diploid is the progeny of two haploids, the two rows that
#'  represent that diploid are copies of the two parental haploid rows
#'  2. When a haploid is the progeny of a diploid, the row that represents that
#'  haploid is the average of the two rows that represent the diploid parent
#'  3. All diagonal elements are equal to 1
#'  4. The matrix is symmetrical
#'
#' @param pedColumns A data.frame with three columns. The first column
#'  has to be the row number and sire and dam columns refer directly to rows
#'  of the pedigree. Parents of founders need to be set to 0 (ZERO). The row
#'  of a child has to be after (i.e. a higher row number) that of its parents.
#'  Indicate a haploid progeny by setting the second column to NA. If both
#'  columns have numbers in them, the progeny is assumed to be diploid.
#'  If an individual has one known and one unknown parent, set the unknown
#'  parent to 0.
#'
#' @return A matrix, \code{hccMat}, the haploid coefficient of coancestry matrix
#'
calcCCmatrixHaploid <- function(pedColumns){
    calcCCmatRow <- function(pedRec){ # Function to process one row of pedColumns
        prog <- pedRec[1]
        sire <- pedRec[2]
        dam <- pedRec[3]
        progM1 <- prog - 1
        if (sire != 0){
            sireRow <- ccMat[sire, 1:progM1] # Non-founder
        } else{
            sireRow <- integer(progM1) # Founder
        }
        if (is.na(dam)){ # One haploid in an SP that is a copy of a GP parent
            ccMat[prog, 1:progM1] <<- ccMat[1:progM1, prog] <<- sireRow
        } else{ # The haploid of a GP produced by a diploid SP
            if (dam != 0){
                damRow <- ccMat[dam, 1:progM1] # Non-founder
            } else{
                damRow <- integer(progM1) # Founder
            }
            ccMat[prog, 1:progM1] <<- ccMat[1:progM1, prog] <<- (sireRow + damRow) / 2
        }
    }#END calcCCmatRow
    
    # calculate Coef Coan matrix here
    hapPed <- makeHaploidPed(pedColumns)
    hapOnePointer <- hapPed$hapOnePointer
    dipTwoPointers <- hapPed$dipTwoPointers
    hapPed <- hapPed$hapPed
    nHap <- nrow(hapPed)
    ccMat <- diag(nHap)
    tmp <- apply(hapPed[2:nHap,], 1, calcCCmatRow)
    return(list(hccMat=ccMat, mccMat=condenseHapCCmat(ccMat, hapOnePointer, dipTwoPointers), hapPed=hapPed, hapOnePointer=hapOnePointer, dipTwoPointers=dipTwoPointers))
}

#' Calculate a haploid pedigree from a mixed pedigree
#'
#'  In the haploid pedigree a diploid is represented by two rows, one for each
#'  of the maternal and paternal haploid parents. The pedigree has three
#'  columns. The first column is the id. It is also the row number. The second
#'  and third columns are pointers to parental haploids.  Each haploid of a
#'  diploid points to only one haploid parent.  A haploid has two pointers, one
#'  to each haploid of its diploid parent.
#'
#' @param pedColumns A data.frame with three columns. The first column
#'  has to be the row number and sire and dam columns refer directly to rows
#'  of the pedigree. Parents of founders need to be set to 0 (ZERO). The row
#'  of a child has to be after (i.e. a higher row number) that of its parents.
#'  Indicate a haploid progeny by setting the second column to NA. If both
#'  columns have numbers in them, the progeny is assumed to be diploid.
#'  If an individual has one known and one unknown parent, set the unknown
#'  parent to 0.
#'
#' @return A matrix, \code{hccMat}, the haploid coefficient of coancestry matrix
#'
makeHaploidPed <- function(pedColumns){
    # The haplotype pedigree to be constructed
    hapPed <- NULL
    
    # What row we are in as we construct hapPed
    runningRowCount <- 1
    
    # Each haploid in pedColumns has one row in hapPed. This points to it.
    hapOnePointer <- NULL
    
    # Each diploid in pedColumns has two rows in hapPed. This points to them.
    dipTwoPointers <- NULL
    makeHapPedRow <- function(pedRow,
    runningRowCount,
    hapOnePointer,
    dipTwoPointers){
        if (is.na(pedRow[3])){
            # It's a haploid
            hapOnePointer <- rbind(hapOnePointer, c(pedRow[1], runningRowCount))
            if (pedRow[2] == 0) return(list(matrix(c(runningRowCount, 0, 0), nrow=1),
            runningRowCount+1,
            hapOnePointer,
            dipTwoPointers))
            dip <- which(dipTwoPointers[,1] == pedRow[2])
            if (length(dip) == 0)
            stop(paste("Diploid parent of haploid", pedRow[1], "not in pedigree"))
            return(list(matrix(c(runningRowCount, dipTwoPointers[dip, 2:3]), nrow=1),
            runningRowCount+1,
            hapOnePointer,
            dipTwoPointers))
        } else{
            # It's a diploid
            dipTwoPointers <- rbind(dipTwoPointers, c(pedRow[1], runningRowCount + 0:1))
            if (pedRow[2] == 0){
                row1 <- c(runningRowCount, 0, NA)
            } else{
                hap <- which(hapOnePointer[,1] == pedRow[2])
                if (length(hap) == 0)
                stop(paste("Haploid parent of diploid", pedRow[1], "not in pedigree"))
                row1 <- c(runningRowCount, hapOnePointer[hap, 2], NA)
            }
            if (pedRow[3] == 0){
                row2 <- c(runningRowCount+1, 0, NA)
            } else{
                hap <- which(hapOnePointer[,1] == pedRow[3])
                if (length(hap) == 0)
                stop(paste("Haploid parent of diploid", pedRow[1], "not in pedigree"))
                row2 <- c(runningRowCount+1, hapOnePointer[hap, 2], NA)
            }
            return(list(rbind(row1, row2),
            runningRowCount+2,
            hapOnePointer,
            dipTwoPointers))
        }
    }
    for (i in 1:nrow(pedColumns)){
        hapPedRow <- makeHapPedRow(pedColumns[i,], runningRowCount, hapOnePointer, dipTwoPointers)
        hapPed <- rbind(hapPed, hapPedRow[[1]])
        runningRowCount <- hapPedRow[[2]]
        hapOnePointer <- hapPedRow[[3]]
        dipTwoPointers <- hapPedRow[[4]]
    }
    rownames(hapPed) <- NULL
    return(list(hapPed=hapPed, hapOnePointer=hapOnePointer, dipTwoPointers=dipTwoPointers))
}

#' Calculate a mixed haploid/diploid coefficient of coancestry matrix from a
#' haploid coefficient of coancestry matrix
#'
#'  In the haploid matrix a diploid is represented by two rows and columns,
#'  one for each of the maternal and paternal gametes. A haploid is represented
#'  by just one row and column
#'
#' @param hccMat The haploid matrix
#' @param hapOnePointer Each haploid has 1 row in hccMat. This points to it.
#' @param dipTwoPointers Each diploid has 2 rows in hccMat. This points to them.
#'
#' @return A matrix, \code{mccMat}, the mixed haploid/diploid
#' coefficient of coancestry matrix
#'
condenseHapCCmat <- function(hccMat, hapOnePointer, dipTwoPointers){
    mccSize <- max(c(hapOnePointer[,1], dipTwoPointers[,1]))
    tmpMat <- matrix(0, nrow=mccSize, ncol=ncol(hccMat))
    mccMat <- matrix(0, nrow=mccSize, ncol=mccSize)
    copyHapToTmp <- function(hapPoint){
        tmpMat[hapPoint[1],] <<- hccMat[hapPoint[2],]
    }
    copyDipToTmp <- function(dipPoint){
        tmpMat[dipPoint[1],] <<- (hccMat[dipPoint[2],] + hccMat[dipPoint[3],])/2
    }
    copyHapToMCC <- function(hapPoint){
        mccMat[,hapPoint[1]] <<- tmpMat[,hapPoint[2]]
    }
    copyDipToMCC <- function(dipPoint){
        mccMat[,dipPoint[1]] <<- (tmpMat[,dipPoint[2]] + tmpMat[,dipPoint[3]])/2
    }
    tmp <- apply(hapOnePointer, 1, copyHapToTmp)
    tmp <- apply(dipTwoPointers, 1, copyDipToTmp)
    tmp <- apply(hapOnePointer, 1, copyHapToMCC)
    tmp <- apply(dipTwoPointers, 1, copyDipToMCC)
    return(mccMat)
}
