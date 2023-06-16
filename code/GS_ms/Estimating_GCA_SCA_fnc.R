# 1. Function calculating GCA
# 2. Cal SCA
# 3. Function running BGLR GCA+SCA model
# 4. Plot out the variance plots
## Two functions, 1 is without any Chks; 2. is if you have Chks
library(stringr)

### 1. If you do not have Checks in the GS model. Y does not have any check rows
CalGCA_noChk  <-  function(Y=Y, colIDy=3, mm.file=mm.file,
                         savefiledir="Onetime1920/Yr19_20/GP1/"){
  print(nrow(Y))
  colIDy <- colIDy

  Y_UniX <- droplevels(Y[!duplicated(Y$Crosses),])  ## Rm duplicated crosses from Y, Linking to each Crosses
  RowNames <- Y$Crosses[!duplicated(Y$Crosses)]

  Y <- Y_UniX  ###  ----> This is for linking to the Individual crosses
  ## Y <- Y  ## If, ----> Linking to full obs. Used for EVD in BGLR

  IDs <- Y[,colIDy]    # the Unique List of P1 OR P2 for ES plots. msZ0 will be based on the length of IDs (aka. P1 or P2)
  print(length(IDs))

  G0  <-  read.table(mm.file,sep=',',header=TRUE,stringsAsFactors=FALSE) # GRM, the first col is names
  name  <-  G0[,1]
  G  <-  G0[,-1]
  colnames(G)  <-  name
  rownames(G)  <-  name

  G <- as.matrix(G)

  Index <- rownames(G)%in%IDs  ##P1:106, Adding index to ensure G has the list of individuals in the IDs (P1 list)
  G <- G[Index,Index]      # subset the G to match up list of P1 names

  phenoNamesFact <- factor(IDs,levels=rownames(G))
    print(length(phenoNamesFact));print(nrow(Y))
  msZ0 <- model.matrix(~-1+phenoNamesFact,data=Y)

    identical(str_split_fixed(colnames(msZ0),"phenoNamesFact",2)[,2],rownames(G))  ## Must be TRUE!!!
  Z <- msZ0   # order of Z individuals is the same as that of G, and plots/Crosses order the same as Y
  colnames(Z) <- stringr::str_remove(colnames(Z),"phenoNamesFact")
      dim(Z)

  GCA <- tcrossprod(tcrossprod(Z,G),Z)  # ZGZ'   oxn,nxn,nxo -> oxo
  rownames(GCA) <- colnames(GCA) <- RowNames

  save(G,GCA,file=paste0(savefiledir,"G.rda"))  ###!!! GP1/
  EVD <- eigen(GCA)
  rownames(EVD$vectors) <- rownames(GCA)
  save(EVD,file=paste0(savefiledir, "EVD.rda")) ###!!! GP1/
  print(dim(msZ0))  # The rownames (row orders) should correspond to those in Y
  #print(dim(chkSpMat))    # 27 * 106 or 27 *121
  print(identical(colnames(Z),rownames(G)))
  print(identical(rownames(Z),rownames(Y)))
  print(dim(GCA))

}


### 2.
CalSCA <- function(G1.file=G1.file, G2.file=G2.file,
                   savefileDir="Onetime1920/Yr19_20/GP1P2/"){
  load(G1.file)
  G1  <-  GCA
  rm(GCA)
  load(G2.file)
  G2  <-  GCA
  GI  <-  G1*G2
  EVD  <-  eigen(GI)
  save(GI, file=paste0(savefileDir,"G.rda"))
  save(EVD, file=paste0(savefileDir,"EVD.rda"))
}


### 3. If you have Checks in the GS model, modify the Z part by adding diagonal matrix
CalGCA <- function(Y=Y, mm.file=mm.file ,colIDy=3, colNam="P1",
                   savefiledir="Onetime1920/Yr19_20/GP1/"){

  colIDy <- colIDy
  colnames(Y)[colIDy] <- colNam      # P1 or P2

  nrowchk <- nrow(Y[Y$crossID=="Check",])      #These are the rows of checks: 33
  dataNHpi_RMchk <- droplevels(Y[!Y$crossID=="Check",]) # Subset ES plots only
  print(dim(dataNHpi_RMchk))
  IDs <- dataNHpi_RMchk[,colIDy]     # the List of P1/P2 for ES plots
  print(head(IDs))

  G0  <-  read.table(mm.file,sep=',',header=TRUE,stringsAsFactors=FALSE)
  name  <-  G0[,1]
  G  <-  G0[,-1]
  colnames(G)  <-  name
  rownames(G)  <-  name
  #  G  <-  2*as.matrix(G) ###### ?????? Diagonal becomes 4!
  G <- as.matrix(G)
  #n <-  dim(G)[1]

  Index <- rownames(G)%in%IDs  ##P1:106, Adding index to ensure G has the list of individuals in the IDs (P1 list)
  G <- G[Index,Index]      # subset the G to match up list of P1 names
  phenoNamesFact <- factor(IDs,levels=rownames(G))   #rownames are sorted alphabetically for the G?
  msZ0 <- model.matrix(~-1 +phenoNamesFact,data=dataNHpi_RMchk)
  ### Z linking P1 and P2 to Crosses, instead of the obs
  chkSpMat <- matrix(0, nrowchk, nrow(G))  # checks are 0s in the Z matrix
  rownames(chkSpMat) <- rownames(Y[Y$crossID=="Check",])
  msZ <- rbind(msZ0,chkSpMat)
  Z <- msZ
  # order of Z individuals is the same as that of G, and plots order the same as Y
  colnames(Z) <- stringr::str_remove(colnames(Z),"phenoNamesFact")

  GCA <- tcrossprod(tcrossprod(Z,G),Z)  # ZGZ' -> oxo / nxn

  save(GCA,file=paste0(savefiledir,"G.rda"))  ###!!! GP1/
  EVD <- eigen(GCA)
  rownames(EVD$vectors) <- rownames(GCA)
  save(EVD,file=paste0(savefiledir,"EVD.rda")) ###!!! GP1/
  print(dim(msZ))  # The rownames (row orders) should correspond to those in Y
  print(dim(chkSpMat))    # 27 * 106 or 27 *121
  print(identical(colnames(Z),rownames(G)))
  print(identical(rownames(Z),rownames(Y)))
  print(dim(GCA))
}


### 3.1 Model_noFMLoc
RunBGLR <- function(YearEffects=TRUE, nIter=50000, burnIn=40000,
                    y=y, testing=testing,
                  Inputfiledir=c("OneTime1920/GP1/",
                                 "OneTime1920/GP2/",
                                 "OneTime1920/GP1P2/"),
                  Outputfiledir){
  library(BGLR)
  nIter   <-  nIter
  burnIn  <-  burnIn

  load(paste0(Inputfiledir[1],"EVD.rda"))
  EVD1 <- EVD
  rm(EVD)
  load(paste0(Inputfiledir[2],"EVD.rda"))
  EVD2 <- EVD
  rm(EVD)
  load(paste0(Inputfiledir[3],"EVD.rda"))
  EVD3 <- EVD
  rm(EVD)


  if (YearEffects){
    ETA <- list(list(~factor(popChk)+factor(Year),data=Y,model="FIXED"),
              list(~factor(line)+factor(block),data=Y,model="BRR"),

              list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
              list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
              list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
    )
  }else{
    ETA <- list(list(~factor(popChk),data=Y,model="FIXED"),
              list(~factor(line)+factor(block),data=Y,model="BRR"),

              list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
              list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
              list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
    )
  }

  fm=BGLR(y=y,
          ETA=ETA,
          nIter=nIter,
          burnIn=burnIn,
          saveAt=paste0(Outputfiledir,"Model_noFMLoc"),
          verbose=TRUE)
  yHat <- fm$yHat
    save(fm,file=paste0(Outputfiledir,"fm.rda"))
  y <- fm$y
  return(list(fm=fm,yHat=fm$yHat))
}

# 3.2 Cor with No FMLoc
predict <- function(testing=testing,yBLUE=yBLUE,Y=Y,fmfiledir){
  load(paste0(fmfiledir,"fm.rda"))
  fm <- fm
  yPred <- fm$ETA[[5]]$u+fm$ETA[[4]]$u+fm$ETA[[3]]$u  #SCA+GCA2+GCA1
  predict <- data.frame(testing,
                      Crosses=Y[c(testing),]$Crosses,
                      yBLUE=yBLUE[testing],
                      yPred=yPred[testing],
                      popChk=Y[c(testing),"popChk"],
                      Year=Y[c(testing),]$Year)  ### !!! gid[testing] is WRONG! # some Crosses are dup between years
  predict <- droplevels(predict)
  write.csv(predict,paste0(fmfiledir,"TP_predict_PP_noFMloc.csv"))
  #predictUniCross <- predict[!duplicated(predict$Crosses),]  #Get the unique Crosses, which may be from 2 years. So no.
  cor <- cor(predict[predict$popChk=="ES",]$yPred,predict[predict$popChk=="ES",]$yBLUE,use="complete")
  return(cor)
}

#### 4.
yHatVarMean_noloc <- function(filedir,vB=2,vGCA1=3,vGCA2=4,vSCA=5,filename="Model_noFMLoc"){
  load(paste0(filedir,"fm.rda"))
  fm <- fm
  varfm <- c(fm$varE,fm$ETA[[vB]]$varB,fm$ETA[[vGCA1]]$varU,fm$ETA[[vGCA2]]$varU,fm$ETA[[vSCA]]$varU)

  varE <- scan(paste0(filedir,filename,"varE.dat"))
  varB <- scan(paste0(filedir,filename,"ETA_",vB,"_varB.dat"))
  varU1 <- scan(paste0(filedir,filename,"ETA_",vGCA1,"_varU.dat"))
  varU2 <- scan(paste0(filedir,filename,"ETA_",vGCA2,"_varU.dat"))
  varU3 <- scan(paste0(filedir,filename,"ETA_",vSCA,"_varU.dat"))

  # plot(varU3,type='o',col=2,cex=.5)
  #### Save the plot!!!!!
  pdf(paste0(filedir,"varE_varGCA1_2_SCA_noLoc_test.pdf"))
  par(mfrow=c(2,2))
  plot(varE,type='o',col=2,cex=.5)
  plot(varU1,type='o',col=1,cex=.5)
  plot(varU2,type='o',col=1,cex=.5)
  plot(varU3,type='o',col=1,cex=.5)
    dev.off()

  varMean <- c(mean(varE),mean(varB),mean(varU1),mean(varU2),mean(varU3))
  names(varMean) <- c("varE","varB_ETA3","varGCA1","varGCA2","varSCA")

  varMean <- rbind(varMean,varfm)
    write.csv(varMean,paste0(filedir,"varMean_noLoc.csv"))
  return(list(varMean=varMean))
}
