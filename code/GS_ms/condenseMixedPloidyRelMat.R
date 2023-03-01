
#' Calculate a mixed haploid/diploid additive relationship matrix from a
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
#' @return A matrix, \code{mRelMat}, the mixed haploid/diploid
#' additive relationship matrix
#'
condenseMixedPloidyRelMat <- function(hccMat, hapOnePointer, dipTwoPointers){
  mrmSize <- max(c(hapOnePointer[,1], dipTwoPointers[,1]))
  tmpMat <- matrix(0, nrow=mrmSize, ncol=ncol(hccMat))
  mrmMat <- matrix(0, nrow=mrmSize, ncol=mrmSize)
  copyHapToTmp <- function(hapPoint){
    tmpMat[hapPoint[1],] <<- hccMat[hapPoint[2],]
  }
  copyDipToTmp <- function(dipPoint){
    tmpMat[dipPoint[1],] <<- hccMat[dipPoint[2],] + hccMat[dipPoint[3],]
  }
  copyHapToMRM <- function(hapPoint){
    mrmMat[,hapPoint[1]] <<- tmpMat[,hapPoint[2]]
  }
  copyDipToMRM <- function(dipPoint){
    mrmMat[,dipPoint[1]] <<- tmpMat[,dipPoint[2]] + tmpMat[,dipPoint[3]]
  }
  tmp <- apply(hapOnePointer, 1, copyHapToTmp)
  tmp <- apply(dipTwoPointers, 1, copyDipToTmp)
  tmp <- apply(hapOnePointer, 1, copyHapToMRM)
  tmp <- apply(dipTwoPointers, 1, copyDipToMRM)
  return(mrmMat)
}
