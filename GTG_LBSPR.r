
DropBox <- "E:/Dropbox"
OutPath <- paste0(DropBox, "/Projects/PackardProject2014/LBSPR_GTGPaper/")
WD <- "E:/GitRepos/LBSPR_ManuscriptCode"

#########################################
# LBSPR  Model - Hordyk et al 2015 ICES #
#########################################

LBSPRSim <- function(StockPars, FleetPars, SizeBins=NULL, P=0.001, Nage=201) {

  MK <- StockPars$MK 
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Beta <- StockPars$FecB 
  MaxSD <- StockPars$MaxSD
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) SizeBins$ToSize <- Linf + MaxSD * SDLinf
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  FM <- FleetPars$FM 
  SL50 <- FleetPars$SL50 
  SL95 <- FleetPars$SL95 
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  x <- seq(from=0, to=1, length.out=Nage) # relative age vector
  EL <- (1-P^(x/MK)) * Linf # length at relative age 
  rLens <- EL/Linf # relative length 
  SDL <- EL * CVLinf # standard deviation of length-at-age
  
  Nlen <- length(LenMids) 
  Prob <- matrix(NA, nrow=Nage, ncol=Nlen)
  Prob[,1] <- pnorm((LenBins[2] - EL)/SDL, 0, 1) # probablility of length-at-age
  for (i in 2:(Nlen-1)) {
    Prob[,i] <- pnorm((LenBins[i+1] - EL)/SDL, 0, 1) - 
		pnorm((LenBins[i] - EL)/SDL, 0, 1)
  }
  Prob[,Nlen] <- 1 - pnorm((LenBins[Nlen] - EL)/SDL, 0, 1)
  
  # Truncate normal dist at MaxSD 
  mat <- array(1, dim=dim(Prob))
  for (X in 1:Nage) {
    ind <- which(abs((LenMids - EL[X]) /SDL[X]) >= MaxSD)
    mat[X,ind] <- 0
  }
  
  Prob <- Prob * mat

  SL <- 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity at length
  Sx <- apply(t(Prob) * SL, 2, sum) # Selectivity at relative age 
  MSX <- cumsum(Sx) / seq_along(Sx) # Mean cumulative selectivity for each age 
  Ns <- (1-rLens)^(MK+(MK*FM)*MSX) # number at relative age in population
  
  Cx <- t(t(Prob) * SL) # Conditional catch length-at-age probablilities  
  Nc <- apply(Ns * Cx, 2, sum) # 
  Pop <- apply(Ns * Prob, 2, sum)
  
  Ml <- 1/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity at length
  Ma <-  apply(t(Prob) * Ml, 2, sum) # Maturity at relative age 
  
  N0 <- (1-rLens)^MK # Unfished numbers-at-age 
  SPR <- sum(Ma * Ns * rLens^Beta)/sum(Ma * N0 * rLens^Beta)
  
  Output <- NULL 
  Output$SPR <- SPR 
  Output$LenMids <- LenMids
  Output$PropLen <- Nc/sum(Nc)
  Output$Pop <- Pop
  
  Output$LCatchFished <- Nc/sum(Nc)
  Output$LPopFished <- Pop
  Output$LCatchUnfished <- apply(N0 * Cx, 2, sum)
  return(Output)
}  

################################################
# GTG LBSPR - Hordyk et al 2016 new manuscript #
################################################

GTGLBSPRSim <- function(StockPars, FleetPars, SizeBins=NULL)  {

  # Assign Variables 
  NGTG <- StockPars$NGTG 
  GTGLinfBy <- StockPars$GTGLinfBy 
  if (!exists("GTGLinfBy")) GTGLinfBy <- NA
  if (is.null(GTGLinfBy)) GTGLinfBy <- NA
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  MaxSD <- StockPars$MaxSD 
  MK <- StockPars$MK 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Walpha <- StockPars$Walpha 
  Wbeta <- StockPars$Wbeta 
  FecB <- StockPars$FecB 
  Steepness <- StockPars$Steepness 
  Mpow <- StockPars$Mpow
  R0 <- StockPars$R0 
 	
  SL50 <- FleetPars$SL50
  SL95 <- FleetPars$SL95 
  MLLKnife <- NA #FleetPars$MLLKnife
  FM <- FleetPars$FM 
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
	
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize

  # Error Catches #
  if (!(exists("NGTG") | exists("GTGLinfBy"))) stop("NGTG or GTGLinfBy must be specified")
  if (!exists("R0")) R0 <- 1E6
  if (is.null(R0)) R0 <- 1E6
  
  # Set up Linfs for the different GTGs
  if (exists("NGTG") & !exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
	GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
  } else  if (!exists("NGTG") & exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
	NGTG <- length(DiffLinfs)
  } else if (exists("NGTG") & exists("GTGLinfBy")) {
    if (!is.na(GTGLinfBy)) {
	  DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
	  NGTG <- length(DiffLinfs)
	} 
	if (is.na(GTGLinfBy)) {
	  DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
	  GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
	}  
  } 
  # Distribute Recruits across GTGS 
  RecProbs <- dnorm(DiffLinfs, Linf, sd=SDLinf) / 
	sum(dnorm(DiffLinfs, Linf, sd=SDLinf)) 
  
  # Length Bins 
  if (is.null(ToSize)) ToSize <- max(DiffLinfs, Linf + MaxSD * SDLinf)
  LenBins <- seq(from=0, by=Linc, to=ToSize)
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=(length(LenBins)-1))

  Weight <- Walpha * LenMids^Wbeta
  
  # Maturity and Fecundity for each GTG 
  L50GTG <- L50/Linf * DiffLinfs # Maturity at same relative size
  L95GTG <- L95/Linf * DiffLinfs # Assumes maturity age-dependant 
  DeltaGTG <- L95GTG - L50GTG
  MatLenGTG <- sapply(seq_along(DiffLinfs), function (X) 
	1.0/(1+exp(-log(19)*(LenMids-L50GTG[X])/DeltaGTG[X])))
  FecLenGTG <- MatLenGTG * LenMids^FecB # Fecundity across GTGs 
  
  VulLen <- 1.0/(1+exp(-log(19)*(LenBins-(SL50+0.5*Linc))/ 
	((SL95+0.5*Linc)-(SL50+0.5*Linc)))) # Selectivity-at-Length
  if (!is.na(MLLKnife)) { # Knife-edge selectivity
    VulLen[LenBins <= MLLKnife] <- 0
	VulLen[LenBins > MLLKnife] <- 1
	SL95 <- SL50 <- NA 
  }

  # Add dome-shaped selectivity curve 
  # Add F-mortality below MLL
  SelLen <- VulLen # Selectivity is equal to vulnerability currently
  
  # Life-History Ratios 
  MKL <- MK * (Linf/(LenBins+0.5*Linc))^Mpow # M/K ratio for each length class
  # Matrix of MK for each GTG
  # MKMat <- sapply(seq_along(DiffLinfs), function(X) 
	# MKL + Mslope*(DiffLinfs[X] - CentLinf))
  MKMat <- matrix(rep(MKL, NGTG), nrow=length(MKL), byrow=FALSE)

  FK <- FM * MK # F/K ratio 
  FKL <- FK * SelLen # F/K ratio for each length class   
  # FkL[Legal == 0] <- FkL[Legal == 0] * DiscardMortFrac 
  ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG
    
  # Set Up Empty Matrices 
  # number-per-recruit at length
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) 
  
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- 
	NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), 
	ncol=NGTG) # number per GTG in each length class 
  # Distribute Recruits into first length class
  NPRFished[1, ] <- NPRUnfished[1, ] <- RecProbs * R0 
  for (L in 2:length(LenBins)) { # Calc number at each size class
    NPRUnfished[L, ] <- NPRUnfished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^MKMat[L-1, ]
    NPRFished[L, ] <- NPRFished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^ZKLMat[L-1, ]
	ind <- DiffLinfs  < LenBins[L]
	NPRFished[L, ind] <- 0
	NPRUnfished[L, ind] <- 0
  } 
  NPRUnfished[is.nan(NPRUnfished)] <- 0
  NPRFished[is.nan(NPRFished)] <- 0
  NPRUnfished[NPRUnfished < 0] <- 0
  NPRFished[NPRFished < 0] <- 0
  
  for (L in 1:length(LenMids)) { # integrate over time in each size class
    NatLUnFishedPop[L, ] <- (NPRUnfished[L,] - NPRUnfished[L+1,])/MKMat[L, ]
    NatLFishedPop[L, ] <- (NPRFished[L,] - NPRFished[L+1,])/ZKLMat[L, ]  
	FecGTGUnfished[L, ] <- NatLUnFishedPop[L, ] * FecLenGTG[L, ]
  }
 
  VulLen2 <- 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity-at-Length
  # print(LenMids)
  # print(c(SL50, SL95))
  # plot(LenMids, VulLen2)
  
  if (!is.na(MLLKnife))  { # Knife-edge selectivity
    VulLen2[LenMids <= MLLKnife] <- 0
	VulLen2[LenMids > MLLKnife] <- 1
	SL95 <- SL50 <- NA 
  }
  
  # points(LenMids, VulLen2, col="red")
  
  # print(cbind(LenMids, VulLen2))
  NatLUnFishedCatch <- NatLUnFishedPop * VulLen2 # Unfished Vul Pop
  NatLFishedCatch <- NatLFishedPop * VulLen2 # Catch Vul Pop
  
  # plot(LenMids, apply(NatLFishedCatch, 1, sum), type="p")
  # matplot(LenMids, (NatLFishedCatch), type="l")
  
  # Expected Length Structure - standardised 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum))
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum))
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))
  
  # Calc SPR
  EPR0 <- sum(NatLUnFishedPop * FecLenGTG) # Eggs-per-recruit Unfished
  EPRf <- sum(NatLFishedPop * FecLenGTG) # Eggs-per-recruit Fished
  SPR <- EPRf/EPR0 
  
  # Equilibrium Relative Recruitment
  recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
  # RelRec/R0 - relative recruitment 
  YPR <- sum(NatLFishedPop  * Weight * VulLen2) * FM 
  Yield <- YPR * RelRec
    
  # Calc Unfished Fitness 
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per Group
  FitPR <- Fit/RecProbs # Fitness per-recruit
  FitPR <- FitPR/median(FitPR)
  ## Debugging
  # plot(FitPR, ylim=c(0,2)) # Should be relatively flat for equal fitness across GTG
   
  # Mslope ignored in this version 
  ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting Mslope 
  Pen <- 0; if (min(MKMat) <= 0 ) Pen <- (1/abs(min(MKMat)))^2 * 1E12 # Penalty for optimising Mslope   
  ObjFun <- ObjFun + Pen
  # print(cbind(Mslope, ObjFun, Pen))

  # Calculate spawning-per-recruit at each size class
  SPRatsize <- cumsum(rowSums(NatLUnFishedPop * FecLenGTG))
  SPRatsize <- SPRatsize/max(SPRatsize)

  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$YPR <- YPR
  Output$LCatchFished <- ExpectedLenCatchFished
  Output$LPopFished <- ExpectedLenPopFished
  Output$LCatchUnfished <- ExpectedLenCatchUnfished
  Output$LPopUnfished <- ExpectedLenPopUnfished
  Output$NatLPopFished <- NatLFishedPop
  Output$NatLPopUnFish <- NatLUnFishedPop
  Output$NatLCatchUnFish <- NatLUnFishedCatch
  Output$NatLCatchFish <- NatLFishedCatch
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- DiffLinfs[2] - DiffLinfs[1]
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$Winf <- Walpha * Linf^Wbeta
  Output$FecLen <- FecLenGTG 
  Output$MatLen <- MatLenGTG 
  Output$SelLen <- SelLen
  Output$MKL <- MKL
  Output$MKMat <- MKMat 
  Output$FKL <- FKL 
  Output$ZKLMat <- ZKLMat 
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$FitPR <- FitPR
  Output$Diff <- range(FitPR)[2] - range(FitPR)[1]
  Output$L50GTG <- L50GTG 
  Output$L95GTG <- L95GTG
  Output$SPRatsize <- SPRatsize
  Output$RelRec <- RelRec
  return(Output)
}

##########################
# Optimisation Functions #
##########################

OptFun <- function(tryFleetPars, LenDat, StockPars, SizeBins=NULL, 
	mod=c("GTG", "LBSPR")) {
  Fleet <- NULL
  Fleet$SL50 <- exp(tryFleetPars[1]) * StockPars$Linf
  Fleet$SL95 <- Fleet$SL50  + (exp(tryFleetPars[2]) * StockPars$Linf)
  Fleet$MLLKnife <- NA
  Fleet$FM <- exp(tryFleetPars[3])
  
  if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, Fleet, SizeBins)
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, Fleet, SizeBins)
  
  LenDat <- LenDat + 1E-15 # add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat)
  predProb <- runMod$LCatchFished 
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb))
  
  # add penalty for SL50 
  trySL50 <- exp(tryFleetPars[1])
  PenVal <- NLL
  Pen <- dbeta(trySL50, shape1=5, shape2=0.01) * PenVal
  if (Pen == 0) Pen <- PenVal * trySL50
  
  # plot(xx, dbeta(xx, shape1=5, shape2=0.01) )
  
  NLL <- NLL+Pen 

  return(NLL)
}

DoOpt <- function(StockPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR")) {
  
  SDLinf <- StockPars$CVLinf * StockPars$Linf
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
 
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf # Starting guesses
  sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
  sFM <- 0.5 
  Start <- log(c(sSL50, sDel, sFM))
  opt <- nlminb(Start, OptFun, LenDat=LenDat, StockPars=StockPars, 
	SizeBins=SizeBins, mod=mod, 
	control= list(iter.max=300, eval.max=400, abs.tol=1E-20))
  
  newFleet <- NULL 
  newFleet$FM <- exp(opt$par[3])
  newFleet$SL50 <- exp(opt$par[1]) * StockPars$Linf 
  newFleet$SL95 <- newFleet$SL50 + exp(opt$par[2]) * StockPars$Linf

  if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, newFleet, SizeBins)
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, newFleet, SizeBins)
  
  Out <- NULL 
  Out$Ests <- c(FM=newFleet$FM, SL50=newFleet$SL50, SL95=newFleet$SL95, 
	SPR=runMod$SPR)
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  return(Out)
}


# Make Stock and Fleet Objects 
StockPars <- NULL 
StockPars$MK <- 1.5
StockPars$NGTG <- 17
StockPars$Linf <- 100
StockPars$CVLinf <- 0.1
StockPars$MaxSD <- 2
StockPars$L50 <- 50 
StockPars$L95 <- 55 
StockPars$FecB <- 3 
StockPars$Walpha <- 0.01
StockPars$Wbeta <- 3 
StockPars$Steepness <- 0.9 
StockPars$Mpow <- 0

SizeBins <- NULL 
SizeBins$Linc <- 5

FleetPars <- NULL
FleetPars$SL50 <- 50 
FleetPars$SL95 <- 55
FleetPars$FM <- 0.8

GTGSim <- GTGLBSPRSim(StockPars, FleetPars, SizeBins)
LBSim <- LBSPRSim(StockPars, FleetPars, SizeBins)

# Compare SPR from two models 
GTGSim$SPR
LBSim$SPR 

# Compare Size Structure - should be identical if F is 0 (or close to)
plot(GTGSim$LenMids, GTGSim$LCatchFished, type="l", lwd=2)
lines(LBSim$LenMids, LBSim$LCatchFished, col="blue", lwd=2)

GTGLenDat <- GTGSim$LCatchFished * 100
LBLenDat <- LBSim$LCatchFished * 100

DoOpt(StockPars, LenDat=GTGLenDat, mod="GTG")
DoOpt(StockPars, LenDat=GTGLenDat, mod="LBSPR")

DoOpt(StockPars, LenDat=LBLenDat, mod="GTG")
DoOpt(StockPars, LenDat=LBLenDat, mod="LBSPR")


##############################
# Simulation over range of F #
# Compare both models        #
##############################
MKSeq <- c(0.5, 1.0, 1.5, 2)
FMSeq <- seq(from=0.1, to=3, length.out=20)
LBIn <- GTGIn <- array(NA, dim=c(length(FMSeq), 2, 4))
LBOut <- GTGOut <- array(NA, dim=c(length(FMSeq), 4, 4))

for (Sp in seq_along(MKSeq)) {
  StockPars$MK <- MKSeq[Sp]
  for (X in seq_along(FMSeq)) {
    FleetPars$FM <- FMSeq[X] 
    
    # Simulate GTG 
    GTGSim <- GTGLBSPRSim(StockPars, FleetPars, SizeBins)
    GTGIn[X,1, Sp] <- FleetPars$FM
    GTGIn[X,2, Sp] <- GTGSim$SPR
    
    # Simulate LBSPR 
    LBSim <- LBSPRSim(StockPars, FleetPars, SizeBins)
    LBIn[X,1, Sp] <- FleetPars$FM
    LBIn[X,2, Sp] <- LBSim$SPR
    
    # Estimate with GTG 
    runopt <- DoOpt(StockPars, LenDat=LBSim$LCatchFished * 100, mod="GTG")
    GTGOut[X,1, Sp] <- runopt$Ests[1]
    GTGOut[X,2, Sp] <- runopt$Ests[4]
    runopt <- DoOpt(StockPars, LenDat=GTGSim$LCatchFished * 100, mod="GTG")
    GTGOut[X,3, Sp] <- runopt$Ests[1]
    GTGOut[X,4, Sp] <- runopt$Ests[4]
    
    # Estimate with LBSPR 
    runopt <- DoOpt(StockPars, LenDat=GTGSim$LCatchFished * 100, mod="LBSPR")
    LBOut[X,1, Sp] <- runopt$Ests[1]
    LBOut[X,2, Sp] <- runopt$Ests[4]
    runopt <- DoOpt(StockPars, LenDat=LBSim$LCatchFished * 100, mod="LBSPR")
    LBOut[X,3, Sp] <- runopt$Ests[1]
    LBOut[X,4, Sp] <- runopt$Ests[4]
    
    print(paste0(X, "/", length(FMSeq)))
    flush.console()
  }
}

Xline <- 3
Yline <- 2.5 
XCex <- 1.25
TexCex <- 1.2
jpeg(paste0(OutPath, "/NewFigures/CompareTwoMods.jpg"), width=160, 
	height=160, res=300, units="mm")
par(mfrow=c(2,2), mar=c(3,1,1,1), oma=c(2,3,0,0))
plot(range(FMSeq), range(FMSeq), type="l", lty=3, bty="l", 
	ylim=c(0, max(LBOut)), axes=FALSE, xlab="", ylab="")
for (X in 1:4) lines(GTGIn[,1,X],  LBOut[,1, X], lwd=2, lty=X)
axis(side=1, label=TRUE)
axis(side=2)
mtext(side=1, expression(italic(F/M)), line=Xline, cex=XCex)
mtext(side=2, bquote("Est" ~ italic(F/M)), line=Yline, cex=XCex)
text(0.13, 4, "a)", cex=TexCex, xpd=NA)

plot(range(FMSeq), range(FMSeq), type="l", lty=3, bty="l", 
	ylim=c(0, max(LBOut)), axes=FALSE, xlab="", ylab="")
for (X in 1:4) lines(LBIn[,1,X],  GTGOut[,1, X], lwd=2, lty=X)
axis(side=1, label=TRUE)
axis(side=2, label=FALSE)
mtext(side=1, expression(italic(F/M)), line=Xline, cex=XCex)
text(0.13, 4, "b)", cex=TexCex, xpd=NA)

plot(c(0,1), c(0,1), type="l", lty=3, bty="l", ylim=c(0, 1), axes=FALSE, 
	xlab="", ylab="")
for (X in 1:4) lines(GTGIn[,2,X],  LBOut[,2, X], lwd=2, lty=X)
axis(side=1, label=TRUE)
axis(side=2, label=TRUE)
mtext(side=1, "SPR", line=Xline, cex=XCex)
mtext(side=2, bquote("Est" ~ SPR), line=Yline, cex=XCex)
text(0.01, 1, "c)", cex=TexCex, xpd=NA)

plot(c(0,1), c(0,1), type="l", lty=3, bty="l", ylim=c(0, 1), 
	axes=FALSE, xlab="", ylab="")
for (X in 1:4) lines(LBIn[,2,X],  GTGOut[,2, X], lwd=2, lty=X)
axis(side=1)
axis(side=2, label=FALSE)
mtext(side=1, "SPR", line=Xline, cex=XCex)
text(0.01, 1, "d)", cex=TexCex, xpd=NA)

legend("bottomright", lty=1:4, lwd=2, legend=MKSeq, 
	title=expression(italic(M/K)), bty="n", cex=1.2)

dev.off()



###################
# Empirical Tests #
###################

NData <- NULL 
EmpPath <- paste0(DropBox,"/Projects/PackardProject2014/LBSPR_GTGPaper/EmpiricalData")

EmpData <- read.csv(paste(EmpPath, "EmpData.csv", sep="/"), header=TRUE, 
  stringsAsFactors = FALSE)
Names <- EmpData[,1]

Species=character()
SPRLB=numeric()
SPRGTG=numeric()
FMLB=numeric()
FMGTG=numeric()
SL50LB=numeric()
SL50GTG=numeric()
SL95LB=numeric()
SL95GTG=numeric()
N=numeric()

par(mfrow=c(3,4))
for (SpNum in 1:12) {# <- 4
  StockPars <- NULL 
  StockPars$MK <- as.numeric(EmpData[which(Names == "MK"), SpNum +1])
  StockPars$NGTG <- as.numeric(EmpData[which(Names == "NGTG"), SpNum +1])
  StockPars$Linf <- as.numeric(EmpData[which(Names == "Linf"), SpNum +1])
  StockPars$CVLinf <- as.numeric(EmpData[which(Names == "CVLinf"), SpNum +1])
  StockPars$MaxSD <- as.numeric(EmpData[which(Names == "MaxSD"), SpNum +1])
  StockPars$L50 <- as.numeric(EmpData[which(Names == "L50"), SpNum +1]) 
  StockPars$L95 <- as.numeric(EmpData[which(Names == "L95"), SpNum +1]) 
  StockPars$FecB <- as.numeric(EmpData[which(Names == "FecB"), SpNum +1]) 
  StockPars$Walpha <- as.numeric(EmpData[which(Names == "Walpha"), SpNum +1])
  StockPars$Wbeta <- as.numeric(EmpData[which(Names == "Wbeta"), SpNum +1]) 
  StockPars$Steepness <- 0.9 
  StockPars$Mpow <- 0
  
  SizeBins <- NULL
  SizeBins$Linc <- as.numeric(EmpData[which(Names == "DatLinc"), SpNum +1]) 
  SizeBins$ToSize <- StockPars$Linf * 1.25
  
  File <- EmpData[which(Names == "DataFile"), SpNum +1]
  EmpDat <- read.csv(paste0(EmpPath, "/", File), header=FALSE, 
    stringsAsFactors = FALSE)
    
  NData[SpNum] <- sum(length(unlist(EmpDat)))
  # }
  
  LenBins <- seq(from=0, to=SizeBins$ToSize, by=SizeBins$Linc)
  LenMids <- seq(from=0.5*SizeBins$Linc, by=SizeBins$Linc, 
  	length.out=(length(LenBins)-1))
  LenDat <- as.vector(table(cut(unlist(EmpDat), LenBins)))	
  
  GTGEsts <- DoOpt(StockPars, LenDat=LenDat, SizeBins=SizeBins, mod="GTG")
  LBEsts <- DoOpt(StockPars, LenDat=LenDat, SizeBins=SizeBins, mod="LBSPR")
  
  Species <- append(Species, File)
  SPRLB <- append(SPRLB, LBEsts$Ests[4])
  SPRGTG <- append(SPRGTG, GTGEsts$Ests[4])
  FMLB <- append(FMLB, LBEsts$Ests[1])
  FMGTG <- append(FMGTG, GTGEsts$Ests[1])
  SL50LB <- append(SL50LB, LBEsts$Ests[2])
  SL50GTG <- append(SL50GTG, GTGEsts$Ests[2]) 
  SL95LB <- append(SL95LB, LBEsts$Ests[3])
  SL95GTG <- append(SL95GTG, GTGEsts$Ests[3])
  N <- append(N, sum(length(unlist(EmpDat))))
  
  tt <- barplot(LenDat, names.arg=LenMids)
  lines(tt, GTGEsts$PredLen, lwd=2)
  lines(tt, LBEsts$PredLen, col="blue", lwd=2)


 print(paste0(SpNum, "/12"))
 flush.console()

}


DF <- data.frame(Species, SPRLB, SPRGTG,
				 FMLB, FMGTG, SL50LB,
				 SL50GTG, SL95LB, SL95GTG,
				 N)
				 

DF

ind <- which(N > 400)
DF[ind,]