# Likelihood in R
## Read in true data

IO_epi <- function(dataPath = "."){
  epiDF_input <- read.csv("./epi_sim.csv",header=T)
  return(epiDF_input)
}

IO_genetic_data_i <- function(ID){
  subject_T_name <- paste0("subject_",ID,"_t_nt.txt")
  
  if(!length(count.fields(subject_T_name))){
    return(list(c(-99),c(-99),c(-99)))
  }
  
  subject_T <- read.csv(subject_T_name,header = F)[,1]
  
  subject_DT <- diff(subject_T)
  
  subject_NT_name <- paste0("subject_",ID,"_nt.txt")
  
  if(!length(count.fields(subject_NT_name))){
    return(list(c(-99),c(-99),c(-99)))
  }
  
  subject_NT <- read.csv(subject_NT_name,header = F)
  
  subject_distance <- rowSums(subject_NT[(1:(nrow(subject_NT)-1)),] != subject_NT[(2:(nrow(subject_NT))),])
  
  return(list(subject_T,subject_DT,subject_distance))
}

IO_source <- function(dataPath){
  src_input <- read.table("infected_source.txt")[,1]
  return(src_input)
}

IO_source_R <- function(dataPath){
  src_input <- read.table("infected_source.txt")[,1] + 1
  return(src_input)
}

IO_genetic_distance <- function(dataPath=".",npop){
  IDs <- 0:(npop-1)
  
  geneticData <- lapply(IDs,IO_genetic_data_i)
  
  return(geneticData)
}

genetic_likelihood_i <- function(ID,geneticData,lambda){
  if((geneticData[[ID]][[2]][1]) < 0){
    return(0)
  }
  
  poissonMeans <- geneticData[[ID]][[2]]*lambda
  
  llhood <- dpois(geneticData[[ID]][[3]],poissonMeans,log=T)
  
  return(sum(llhood,na.rm=T))
}

genetic_likelihood <- function(geneticData,paramsDF){
  npop <- length(geneticData)
  return(sapply(1:npop,genetic_likelihood_i,geneticData=geneticData,lambda=paramsDF$lambda))
}

distance_matrix <- function(longVec,latVec){
  return(as.matrix(dist(cbind(longVec,latVec))))
}

exp_kernel <- function(distance,kappa){
  return(exp(-kappa*distance))
}

accumulate_exposure_i <- function(tNow,ID,paramDF,epiDF,distMat,tStart = 0){
  infectTime <- pmin(epiDF[,"t_r"],tNow) - pmax(epiDF[,"t_i"],tStart)
  
  infectInterval <- infectTime
  infectInterval[infectTime <= 0] <- 0
  infectInterval[ID] <- 0
  
  bgInf <- (tNow - tStart)*paramDF[,"alpha"]
  
  secondInf <- infectInterval*exp_kernel(distMat[,ID],paramDF[,"kappa"])*paramDF[,"beta"]
  
  return(sum(bgInf,secondInf))
}

accumulate_exposure_t <- function(tNows,ID,paramDF,epiDF,distMat,tStart = 0){
  return(sapply(tNows,accumulate_exposure_i,ID=ID,paramDF=paramDF,epiDF=epiDF,distMat=distMat,tStart=tStart))
}

instant_exposure_i <- function(ID,src,paramDF,distMat){
  exposure <- 1
  if(src > paramDF[,"n"]){
    exposure <- paramDF[,"alpha"]
  } else if(src == ID){
    exposure <- 1
  } else if(src > 0){
    exposure <- paramDF[,"beta"] * exp_kernel(distMat[ID,src],paramDF[,"kappa"])
  }
  return(log(exposure))
}

instant_exposure <- function(IDs,srcs,paramDF,distMat){
  exposure <- rep(1,paramDF[,"n"])
  secondSrc <- srcs[(srcs > 0) & (srcs <= paramDF[,"n"])]
  secondID <- IDs[(srcs > 0) & (srcs <= paramDF[,"n"])]
  distCall <- cbind(secondID,secondSrc)
  dists <- distMat[distCall]
  exposure[(srcs > 0) & (srcs <= paramDF[,"n"])] <- paramDF[,"beta"]*exp_kernel(dists,paramDF[,"kappa"])
  exposure[srcs > paramDF[,"n"]] <- paramDF[,"alpha"]
  
  return(log(exposure))
}

xU_likelihood <- function(epiDF,paramDF,distMat){
  t_E <- pmin(epiDF[,"t_e"],paramDF[,"tmax"])
  return(rowSums(sapply(1:(paramDF[,"n"]),accumulate_exposure_t,tNows=t_E,paramDF=paramDF,epiDF=epiDF,distMat=distMat)))
}

xE_likelihood <- function(epiDF,srcs,paramDF,distMat){
  instantExposure <- instant_exposure(IDs=1:(paramDF[,"n"]),srcs=srcs,paramDF=paramDF,distMat=distMat)
  return(instantExposure)
}

xEnI_likelihood <- function(epiDF,paramDF){
  which_xEnI <- (epiDF[,"t_i"] > paramDF[,"tmax"]) & (epiDF[,"t_e"] <= paramDF[,"tmax"])
  infect_interval <- which_xEnI*(paramDF[,"tmax"] - epiDF[,"t_e"])
  F_E <- pgamma(infect_interval,shape=paramDF[,"a"],scale=paramDF[,"b"],log.p = T,lower.tail = F)
  F_E[!which_xEnI] <- 0
  
  return(F_E)
}

xI_likelihood <- function(epiDF,paramDF){
  which_xI <- (epiDF[,"t_i"] <= paramDF[,"tmax"]) & (epiDF[,"t_e"] <= paramDF[,"tmax"])
  infect_interval <- which_xI*(epiDF[,"t_i"] - epiDF[,"t_e"])
  f_I <- dgamma(infect_interval,shape=paramDF[,"a"],scale=paramDF[,"b"],log = T)
  f_I[!which_xI] <- 0
  
  return(f_I)
}

xInR_likelihood <- function(epiDF,paramDF){
  which_xInR <- (epiDF[,"t_r"] > paramDF[,"tmax"]) & (epiDF[,"t_i"] <= paramDF[,"tmax"])
  infect_interval <- which_xInR*(paramDF[,"tmax"] - epiDF[,"t_i"])
  F_I <- pweibull(infect_interval,shape=paramDF[,"c"],scale=paramDF[,"d"],log.p = T,lower.tail = F)
  F_I[!which_xInR] <- 0
  
  return(F_I)
}

xR_likelihood <- function(epiDF,paramDF){
  which_xR <- (epiDF[,"t_r"] <= paramDF[,"tmax"]) & (epiDF[,"t_i"] <= paramDF[,"tmax"])
  infect_interval <- which_xR*(epiDF[,"t_r"] - epiDF[,"t_i"])
  f_R <- dgamma(infect_interval,shape=paramDF[,"a"],scale=paramDF[,"b"],log = T)
  f_R[!which_xR] <- 0
  
  return(f_R)
}

test_accumulate_exposure <- function(){
  epidat <- data.frame("coor_x" = c(0,2,3,3,4,5),"coor_y" = c(0,3,9,6,7,1), "t_i" = c(0,1,7,6,50,900000), "t_r" = c(50,5,15,9000000,70,9000000))
  distmat <- distance_matrix(epidat$coor_x,epidat$coor_y)
  paramtest <- data.frame("alpha"=0.0001,"beta"=0.7,"kappa"=1.7,"lambda"=0.2)
  
  accumulate_exposure(10,1:10,paramtest,epidat,distmat)
  
  print(distmat[1,] == sqrt(c(0,2,3,3,4,5)^2 + c(0,3,9,6,7,1)^2))
  
  kerneltest <- exp(-sqrt(c(0,2,3,3,4,5)^2 + c(0,3,9,6,7,1)^2)*1.7)
  
  print(exp_kernel(distmat[1,],paramtest$kappa) == kerneltest)
  
  intervaltest <- pmin(epidat$t_r,10) - epidat$t_i
  intervaltest[intervaltest <= 0] <- 0
  intervaltest[1] <- 0
  
  print(sum(intervaltest*0.7*kerneltest,0.0001*10) == accumulate_exposure(10,1,paramtest,epidat,distmat))
}

lh_square <- function(paramDF,epiDF,geneDF,src){
  distMat <- distance_matrix(epiDF[,"coor_x"],epiDF[,"coor_y"])
  
  lhSquare <- data.frame(l_U=xU_likelihood(epiDF,paramDF,distMat),
                         l_E=xE_likelihood(epiDF,src,paramDF,distMat),
                         l_EnI=xEnI_likelihood(epiDF,paramDF),
                         l_I=xI_likelihood(epiDF,paramDF),
                         l_InR=xInR_likelihood(epiDF,paramDF),
                         l_R=xR_likelihood(epiDF,paramDF),
                         l_gene=genetic_likelihood(geneData,params))
  
  return(lhSquare)
}

likelihood <- function(paramDF,epiDF,geneDF,src){
  return(sum(lh_square(paramDF,epiDF,geneDF,src)))
}

lh_square_dir <- function(paramDF,simOutputPath="."){
  epiDF <- IO_epi(simOutputPath)
  src <- IO_source_R(simOutputPath)
  npop <- nrow(epiDF)
  geneDF <- IO_genetic_distance(dataPath = simOutputPath, npop = npop)
  
  return(lh_square(paramDF=paramDF,epiDF=epiDF,geneDF=geneDF,src=src))
}

likelihood_dir <- function(paramDF,simOutputPath="."){
  epiDF <- IO_epi(simOutputPath)
  src <- IO_source_R(simOutputPath)
  npop <- nrow(epiDF)
  geneDF <- IO_genetic_distance(dataPath = simOutputPath, npop = npop)
  
  return(likelihood(paramDF=paramDF,epiDF=epiDF,geneDF=geneDF,src=src))
}
