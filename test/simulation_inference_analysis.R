#!/usr/bin/env Rscript
library(BORIS2)
library(ggplot2)

seed <- 108
cluster <- T
params_file <- T
long_run <- T
simonly <- F

subsampling_pct <- 1

kernel_location <- F

npop <- 150
nbase <- 8000
tmax <- 104

alphaval <- 0.0001
betaval <- 8
kappaval <- 0.02
aval <- 10
bval <- 0.5
cval <- 2
dval <- 2
mu1val <- 0.002
mu2val <- 0.0005
pval <- 0.0001

maxx <- 2000
maxy <- 2000

maxx <- 2000
maxy <- 2000

kerneltype <- 'exponential'

set.seed(seed=seed)

## Check if the script is running on the cluster
whichdir <- getwd()
if(substr(whichdir,2,6) == "Users"){
  cluster <- F
}


if(cluster & params_file){
  args = commandArgs(trailingOnly=TRUE)
  filename <- args[1]
  fileline <- as.numeric(args[2])
  
  print(paste0("cluster = t, param_file = ",filename," line ",fileline))
  
  param_settings <- read.table(filename,header=T)
  
  if(fileline > nrow(param_settings)){
    # HALT
    UNDEFINED()
  }
  
  param_values <- as.numeric(param_settings[fileline,])
  
  print("Parameter values:")
  print(param_values)
  
  seed <- param_values[1]
  alphaval <- param_values[2]
  betaval <- param_values[3]
  kappaval <- param_values[4]
  aval <- param_values[5]
  bval <- param_values[6]
  cval <- param_values[7]
  dval <- param_values[8]
  mu1val <- param_values[9]
  mu2val <- param_values[10]
  maxx <- param_values[11]
  maxy <- param_values[12]
  tmax <- param_values[13]
  npop <- param_values[14]
  nbase <- param_values[15]
  subsampling_pct <- param_values[16]
  kerneltype <- "exponential"
  
  print(paste0("Setting seed: ",seed))
  set.seed(seed)
  
  seedpaste <- seed
  linepaste <- fileline
  
  if(seed < 10){
    seedpaste <- paste0("0",seed)
  }
  if(fileline < 10){
    linepaste <- paste0("0",fileline)
  }
  
  dirname <- paste0("/projects/lau_projects/phylodynamics_hannah/sim_test_line_",linepaste,"_seed_",seedpaste)
  
  if(!dir.exists(dirname)){
    dir.create(dirname)
  }
  setwd(dirname)
} 

if(!is.na(seed) & !cluster){
  seedpaste <- seed
  
  if(seed < 10){
    seedpaste <- paste0("0",seed)
  }
  
  dirname <- paste0("~/Documents/Research/Phylodynamics/BORIS2_Results/sim_test_seed_",seedpaste)
  if(!dir.exists(dirname)){
    dir.create(dirname)
  }
  setwd(dirname)
}

if(cluster & !params_file){
  args = commandArgs(trailingOnly=TRUE)
  seed <- as.numeric(args[1])
  set.seed(seed)
  
  seedpaste <- seed
  
  if(seed < 10){
    seedpaste <- paste0("0",seed)
  }
  
  dirname <- paste0("/projects/lau_projects/phylodynamics_hannah/sim_test_seed_",seedpaste)
  if(!dir.exists(dirname)){
    dir.create(dirname)
  }
  setwd(dirname)
}

ti_known <- T
te_known <- T
ab_known <- T
cd_known <- T
abk_known <- T
ndiff_known <- T
ndiff_true_init <- T
source_known <- T
lambda_init_scale <- 1

run_params <- data.frame(seed=seed,
                         npop=npop,
                         nbase=nbase,
                         tmax=tmax,
                         ti_known=ti_known,
                         te_known=te_known,
                         ab_known=ab_known,
                         cd_known=cd_known,
                         abk_known=abk_known,
                         ndiff_known=ndiff_known,
                         ndiff_true_init=ndiff_true_init,
                         source_known=source_known,
                         lambda_init_scale=lambda_init_scale,
                         subsampling_pct=subsampling_pct)

data(sim.epi.input)

sim.epi.data <- sim.epi.input

data(sim.moves.input)

sim.epi.data <- data.frame(k=0:(npop-1),
                           coor_x=runif(npop,0,maxx),
                           coor_y=runif(npop,0,maxy),
                           t_e=c(0,rep(9e+06,npop-1)),
                           t_i=rep(9e+06,npop),
                           t_r=rep(9e+06,npop),
                           ftype0=sample(c(0,1),npop,replace=T),
                           herdn=sample(min(sim.epi.input$herdn):max(sim.epi.input$herdn),
                                        npop,
                                        replace=T),
                           status=c(2,rep(1,npop-1)))

sim.epi.data$ftype1=ifelse(sim.epi.data$ftype0==1,
                           0,
                           sample(c(0,1),npop,replace=T))
sim.epi.data$ftype2=ifelse(sim.epi.data$ftype0==0 & sim.epi.data$ftype1==0,
                           1,
                           0)

sim.epi.data <- sim.epi.data[,names(sim.epi.input)]

if(long_run==F){
  niter = 20
  nfrequ = 1
  noutputsource = 1
  noutputgm = 1
  ncout = 1
} else{
  niter = 300000
  nfrequ = max(round(npop/3),1)
  noutputsource = max(round(niter/10000),1)
  noutputgm = max(round(niter/1000),1)
  ncout = max(round(niter/1000),1)
}

nburn = round(niter/10)

# 'alpha' = 7e-4,
# 'beta' = 0.5,
# 'p_ber' = 0.1,
# 'mu_1' = 2e-05,
# 'mu_2' = 2e-06,

## ----sim_setup2--------------------------------------------------------------------------------
para.key <- data.frame('alpha' = alphaval,
                       'beta' = betaval,
                       'mu_1' = mu1val,
                       'mu_2' = mu2val,
                       'a' = aval,
                       'b' = bval,
                       'c' = cval,
                       'd' = dval,
                       'k_1' = kappaval,
                       'p_ber' = pval,
                       'phi_inf1' = 3,
                       'phi_inf2' = 1.5,
                       'rho_susc1' = 0.4,
                       'rho_susc2' = 2,
                       'nu_inf' = 0.2,
                       'tau_susc'= 0.1,
                       'beta_m'= 1.0)

simulate_location_kernel <- function(npop,beta_t,kappa_t,xmax,ymax){
  x_coord <- rep(NA,npop)
  y_coord <- rep(NA,npop)
  
  distances <- rexp(npop,rate=kappa_t)
  angles <- runif(npop,min=0,max=2*pi)
  
  #index
  x_coord[1] <- runif(1,min=xmax/4,max=(3*xmax)/4)
  y_coord[1] <- runif(1,min=ymax/4,max=(3*ymax)/4)
  
  for(i in 2:npop){
    prev_x <- x_coord[(i-1)]
    prev_y <- y_coord[(i-1)]
    
    prev_dist <- distances[i]
    prev_angle <- angles[i]
    
    add_x <- prev_dist*cos(prev_angle)
    add_y <- prev_dist*sin(prev_angle)
    
    new_x <- prev_x + add_x
    new_y <- prev_y + add_y
    
    if((new_x < 0 | new_x > xmax) | (new_y < 0 | new_y > ymax)){
      prev_angle <- -angles[i]
      
      add_x <- prev_dist*cos(prev_angle)
      add_y <- prev_dist*sin(prev_angle)
      
      new_x <- prev_x + add_x
      new_y <- prev_y + add_y
    }
    while((new_x < 0 | new_x > xmax) | (new_y < 0 | new_y > ymax)){
      prev_dist <- rexp(1,rate=kappa_t)
      prev_angle <- runif(1,min=0,max=2*pi)
      
      add_x <- prev_dist*cos(prev_angle)
      add_y <- prev_dist*sin(prev_angle)
      
      new_x <- prev_x + add_x
      new_y <- prev_y + add_y
    }
    
    x_coord[i] <- new_x
    y_coord[i] <- new_y
  }
  
  location <- data.frame(x=x_coord,y=y_coord)
  
  return(location)
}

if(kernel_location == T){
  location_data <- simulate_location_kernel(npop,para.key$beta,para.key$k_1,maxx,maxy)
  
  sim.epi.data$coor_x <- location_data$x
  sim.epi.data$coor_y <- location_data$y
}

## ----sim_setup3--------------------------------------------------------------------------------
pars.aux <- data.frame('n' = nrow(sim.epi.data),
                       'seed' = seed,
                       'n_base' = nbase,
                       'n_seq' = 5,
                       't_max' = tmax,
                       'unassigned_time' =  9e+6,
                       'sample_range' = 10,
                       'partial_seq_out' = 0,
                       'n_base_part' = nbase,
                       'n_index' = 1,
                       'coord_type' = 'cartesian',
                       'kernel_type' = kerneltype,
                       'latent_type' = 'gamma',
                       'opt_k80' = 0,
                       'opt_betaij' = 0,
                       'opt_mov' = 0,
                       'n_mov' = 60,
                       stringsAsFactors = F)

if(!dir.exists("./Sim_Test")){
  dir.create("./Sim_Test")
}

setwd("./Sim_Test")

## ----sim---------------------------------------------------------------------------------------
sim.out<-sim(epi.inputs = sim.epi.data,
             moves.inputs = sim.moves.input,
             parsKey = para.key,
             parsAux = pars.aux,
             inputPath = "./inputs",
             outputPath = "./outputs")

pdf(file="outbreak_map.pdf")
plot(sim.out$epi.sim[,"coor_x"],sim.out$epi.sim[,"coor_y"],col=ifelse(sim.out$epi.sim[,"t_e"]<=tmax,"black","grey"),pch=19)
dev.off()

print(paste0("Outbreak size: ",sum(sim.out$infected_source>=0)))

pdf(file="latent_period.pdf")
hist(sim.out$epi.sim[which(sim.out$epi.sim[,"t_i"] < tmax),"t_i"] - sim.out$epi.sim[which(sim.out$epi.sim[,"t_i"] < tmax),"t_e"],main="Latent period")
hist(sim.out$epi.sim[which(sim.out$epi.sim[,"t_r"] < tmax),"t_r"] - sim.out$epi.sim[which(sim.out$epi.sim[,"t_r"] < tmax),"t_i"],main="Infectious period")
dev.off()

# Define some summary functions
arrow_outbreak_plot <- function(epi.data,true.source.c,t.max=100){
  library(ggplot2)
  library(dplyr)
  
  epi.data.df <- as.data.frame(epi.data)
  true_source_r <- true.source.c + 1
  true_source_r[true_source_r > nrow(epi.data) | true_source_r < 0] <- NA
  epi.data.df$arrow_start_x <- epi.data.df$coor_x[true_source_r]
  epi.data.df$arrow_start_y <- epi.data.df$coor_y[true_source_r]
  
  is.infected <- ifelse(epi.data.df$t_e <= t.max,"yes","no")
  return(epi.data.df %>% ggplot(aes(x=coor_x,y=coor_y,alpha=0.83)) +
           geom_point(mapping=aes(colour=is.infected)) +
           scale_color_manual(values=c("no"="grey85","yes"="white")) +
           geom_point(data=subset(epi.data.df,is.infected=="yes")) +
           theme_classic() +
           geom_point(data=subset(epi.data.df,true_source_r==10000),shape=3) +
           geom_segment(aes(x=arrow_start_x,y=arrow_start_y,xend=coor_x,yend=coor_y),arrow=arrow(length=unit(0.2,'cm')),color="blue") +
           guides(colour="none",is.infected="none",alpha="none"))
}

pressure_dist_plot_infection <- function(kappa,beta,coor_x,coor_y,true_source){
  distmat <- as.matrix(dist(cbind(coor_x,coor_y)))
  
  npop <- length(coor_x)
  
  true_source_r <- true_source + 1
  
  truedist <- rep(NA,length(true_source_r))
  
  for(i in 1:length(true_source_r)){
    if(true_source_r[i] <= npop & true_source_r[i] > 0){
      truedist[i] <- distmat[i,true_source_r[i]]
    }
  }
  
  eq = function(distance){beta*exp(-kappa*distance)}
  
  par(mfrow=c(3,1))
  
  hist(distmat,xlab="Distance, all",xlim=c(0,max(distmat)))
  hist(truedist,xlab="Distance, infection source",xlim=c(0,max(distmat)))
  curve(eq, from=0, to=max(distmat,na.rm=T), xlab="Distance", ylab="Infectious Pressure (B*K)")
  
  par(mfrow=c(1,1))
}

alpha_beta_kernel_comparison <- function(kappa,beta,alpha,coor_x,coor_y,true_source){
  distmat <- as.matrix(dist(cbind(coor_x,coor_y)))
  
  npop <- length(coor_x)
  
  true_source_r <- true_source + 1
  
  truedist <- rep(NA,length(true_source_r))
  
  for(i in 1:length(true_source_r)){
    if(true_source_r[i] <= npop & true_source_r[i] > 0){
      truedist[i] <- distmat[i,true_source_r[i]]
    }
  }
  
  kernel_mat <- beta*exp(-kappa*distmat)
  truekernel <- beta*exp(-kappa*truedist)
  
  eq = function(distance){beta*exp(-kappa*distance)}
  
  par(mfrow=c(3,1))
  
  hist(kernel_mat,xlab="Kernel, all")
  abline(v=alpha,col="blue")
  hist(truekernel,xlab="Kernel, infection source")
  abline(v=alpha,col="blue")
  curve(eq, from=0, to=max(distmat,na.rm=T), xlab="Distance", ylab="Infectious Pressure (B*exp(-k*d))")
  
  par(mfrow=c(1,1))
}

outbreak_plot <- arrow_outbreak_plot(sim.out$epi.sim,sim.out$infected_source,tmax)
pdf(file="outbreak_plot.pdf",width=8,height=8)
print(outbreak_plot)
dev.off()

pdf(file="kernel_summary.pdf",width=7,height=10)
pressure_dist_plot_infection(kappaval,betaval,sim.out$epi.sim[,"coor_x"],sim.out$epi.sim[,"coor_y"],sim.out$infected_source)
alpha_beta_kernel_comparison(kappaval,betaval,alphaval,sim.out$epi.sim[,"coor_x"],sim.out$epi.sim[,"coor_y"],sim.out$infected_source)
dev.off()


## Calculate number infected, see if it's good
print(sum(read.csv("./outputs/infected_source.txt",header=F)[,1] >= 0))
print(sum(read.csv("./outputs/infected_source.txt",header=F)[,1] >= 0)/nrow(sim.epi.data))

## Calculate number infected from bg
print(sum(read.csv("./outputs/infected_source.txt",header=F)[,1] == 9999))
print(sum(read.csv("./outputs/infected_source.txt",header=F)[,1] == 9999)/nrow(sim.epi.data))

# Function to summarize SEIR thru time
countSEIRt <- function(tval,te,ti,tr){
  ns <- sum(te > tval & ti > tval & tr > tval)
  ne <- sum(te <= tval & ti > tval & tr > tval)
  ni <- sum(te <= tval & ti <= tval & tr > tval)
  nr <- sum(te <= tval & ti <= tval & tr <= tval)
  
  return(c(ns,ne,ni,nr))
}

countSEIR <- function(te,ti,tr,tmax.outbreak=100){
  tpoints <- seq(from=0,to=tmax.outbreak,length.out=200)
  ns <- rep(0,length(tpoints))
  ne <- rep(0,length(tpoints))
  ni <- rep(0,length(tpoints))
  nr <- rep(0,length(tpoints))
  
  for(i in 1:length(tpoints)){
    countcat <- countSEIRt(tpoints[i],te,ti,tr)
    
    ns[i] <- countcat[1]
    ne[i] <- countcat[2]
    ni[i] <- countcat[3]
    nr[i] <- countcat[4]
  }
  
  npop <- length(te)
  
  return(data.frame(time=tpoints,NS=ns,NE=ne,NI=ni,NR=nr))
}

sim_data <- read.csv("./outputs/epi_sim.csv")

seir_count <- countSEIR(sim_data$t_e,sim_data$t_i,sim_data$t_r,tmax)

plot(seir_count$time,seir_count$NS,main="SEIR counts",xlab="Time",ylab="N",col="darkgreen",type="l",xlim=c(0,tmax),ylim=c(0,length(sim_data$t_e)))
lines(seir_count$time,seir_count$NE,col="darkblue")
lines(seir_count$time,seir_count$NI,col="darkred")
lines(seir_count$time,seir_count$NR,col="darkgrey")
legend("right",col=c("darkgreen","darkblue","darkred","darkgrey"),lty=1,legend=c("S","E","I","R"))

pdf(file="SEIR_count.pdf",width=8,height=8)
plot(seir_count$time,seir_count$NS,main="SEIR counts",xlab="Time",ylab="N",col="darkgreen",type="l",xlim=c(0,tmax),ylim=c(0,length(sim_data$t_e)))
lines(seir_count$time,seir_count$NE,col="darkblue")
lines(seir_count$time,seir_count$NI,col="darkred")
lines(seir_count$time,seir_count$NR,col="darkgrey")
legend("right",col=c("darkgreen","darkblue","darkred","darkgrey"),lty=1,legend=c("S","E","I","R"))
dev.off()


if(simonly){
  UNDEFINED()
}

dir.create("../Infer_Test")

## ----sim_out1, eval=F--------------------------------------------------------------------------
## sim.out$epi.sim[1:6,]


## ----sim_out2, echo=F--------------------------------------------------------------------------
tmp<-sim.out$epi.sim[1:6,]
tmp[,4:6]<-round(tmp[,4:6],4)
tmp


## ----sim_out3----------------------------------------------------------------------------------
head(sim.out$infected_source, 14)


## ----sim_out4----------------------------------------------------------------------------------
which(sim.out$infected_source == -99)


## ----sim_out5----------------------------------------------------------------------------------
inf <- which(sim.out$infected_source != -99)
length(inf)/length(sim.out$infected_source)


## ----sim_out6----------------------------------------------------------------------------------
sim.out$sampled_perct


## ----sim_out7----------------------------------------------------------------------------------
head(floor(sim.out$t_sample), 20)


## ----sim_out8----------------------------------------------------------------------------------
length(which(sim.out$t_sample[inf] != 9e6))/length(inf)


## ----sim_apeout--------------------------------------------------------------------------------
library(ape)

#read in the sequences from individual k=0
fn <- paste0(sim.out$outputPath, "subject_0_nt.txt")
seqs1 <- as.matrix(read.csv(fn, header=F))

#the sequences for this individul are stored in 19 rows
# nucleotides are coded as 1=A, 2=G, 3=T, 4=C
dim(seqs1)

#convert the sequences to their nucleotide base letters
seqs1[seqs1==1]<-"a"
seqs1[seqs1==2]<-"g"
seqs1[seqs1==3]<-"t"
seqs1[seqs1==4]<-"c"

#name the rows, which become the tip labels in the fasta file.
#for this we will import the timings associated with each sequence
fn <- paste0(sim.out$outputPath, "subject_0_t_nt.txt")
seqs1_t <- as.numeric(unlist(read.table(fn)))

rownames(seqs1) <- paste0("k0_", round(seqs1_t,3))

#write these out to a fasta file
library(ape)
write.dna(seqs1, file='seqs1.fasta', format = 'fasta', nbcol = -1)

#inspect snps in these data
seq.tab <- apply(seqs1, MARGIN = 2, table)
snp.v<-sapply(seq.tab, length)
snps<-as.numeric(which(snp.v>1))
# seqs1[1:8,snps]

#use the epidemiological data from the previous simulation run
epi.data<-as.data.frame(sim.out$epi.sim)

#make sure indexing is correct
epi.data$k<-1:nrow(epi.data)-1

#produce the required initial values for the temporal fields.
#assume that exposure has occurred 24 hours before the onset of infectiousness:
epi.data$t_e <- epi.data$t_i - 1

if(te_known == T){
  epi.data$t_e <- sim.out$epi.sim[,"t_e"]
}

epi.data$t_e[epi.data$t_e > 10000] <- pars.aux$unassigned_time
epi.data$t_r <- round(epi.data$t_r, 0)

#wrangle the farm type variables into the required format
epi.data$ftype <- NA
epi.data$ftype[epi.data$ftype0 == 1] <- 0
epi.data$ftype[epi.data$ftype1 == 1] <- 1
epi.data$ftype[epi.data$ftype2 == 1] <- 2

#keep just the desired columns, re-ordered as required
epi.data<-epi.data[,c('k','coor_x','coor_y','t_e','t_i','t_r','ftype','herdn')]
head(epi.data)


## ----inf_setup2--------------------------------------------------------------------------------
#for each individual propose an initial source for the tree used to initialised the MCMC inference

#this could be assigned manually based on available information and/or opinion. Or, as in the below example randomly sampled from those individuals with earlier onset times than each individual.

epi.data$initial_source<-NA

for (i in 1:nrow(epi.data)){
  
  epi.data$k[i]
  
  #identify the pool of possible sources with earlier initial exposure times
  pool_source<- epi.data$k [which(epi.data$t_e<epi.data$t_e[i] & epi.data$t_r>epi.data$t_e[i])]
  
  #randomly assign an infection source is there are any in the pool
  if (length(pool_source)>1) epi.data$initial_source[i] <- sample(pool_source,1)
  if (length(pool_source)==1) epi.data$initial_source[i] <- pool_source[1]
  
  #otherwise assign 9999 as the initial source, to indicate that it is either external to the dataset or unknown
  if (length(pool_source)<1) epi.data$initial_source[i] <- 9999
  
}
epi.data$initial_source

#use the movement data from previously
moves.inputs<-sim.out$moves.inputs
head(moves.inputs)


## ----inf_setup3--------------------------------------------------------------------------------
pars.aux <- data.frame('n' = nrow(epi.data),
                       'kernel_type' = kerneltype,
                       'coord_type' = 'cartesian',
                       't_max' = tmax,
                       'unassigned_time' =  9e+6,
                       'processes' = 1,
                       'n_seq' = 5,
                       'n_base' = nbase,
                       'n_iterations' = 1e5,
                       'n_frequ' = 10,
                       'n_output_source' = 1000,
                       'n_output_gm' = 2000,
                       'n_cout' = 1000,
                       'opt_latgamma' = 0,
                       'opt_k80' = 0,
                       'opt_betaij' = 0,
                       'opt_ti_update' = 0,
                       'opt_mov' = 0,
                       stringsAsFactors = F)

#once `t_max` is set, right censor all recovery times
epi.data$t_r[epi.data$t_r>pars.aux$t_max]<-pars.aux$t_max

# 'alpha' = 0.1
## ----inf_setup4--------------------------------------------------------------------------------
para.key.inits <- data.frame('alpha' = para.key$alpha,
                             'beta' = para.key$beta*1.5,
                             'lat_mu' = para.key$a,
                             'lat_var' = para.key$b,
                             'c' = para.key$c*2,
                             'd' = para.key$d*2,
                             'k_1' = para.key$k_1*2,
                             'mu_1' = para.key$mu_1*0.5,
                             'mu_2' = para.key$mu_2*0.5,
                             'p_ber' = para.key$p_ber*2,
                             'phi_inf1' = 3,
                             'phi_inf2' = 1.5,
                             'rho_susc1' = 0.4,
                             'rho_susc2' = 2,
                             'nu_inf' = 0.2,
                             'tau_susc'= 0.1,
                             'beta_m'= 1)

# if(ndiff_known & ndiff_true_init){
#   para.key.inits$mu_1 = para.key$mu_1
#   para.key.inits$mu_2 = para.key$mu_2
# }

if(ab_known){
  para.key.inits$lat_mu = para.key$a
  para.key.inits$lat_var = para.key$b^2
}

if(cd_known){
  para.key.inits$c = para.key$c
  para.key.inits$d = para.key$d
}

if(abk_known){
  para.key.inits$alpha = para.key$alpha
  para.key.inits$beta = para.key$beta
  para.key.inits$k_1 = para.key$k_1
}


## ----inf_setup5--------------------------------------------------------------------------------
para.priors <- data.frame('t_range' = 20,
                          't_back' = (para.key$a*para.key$b^2)*2,
                          't_bound_hi' = 40,
                          'rate_exp_prior' = 0.001,
                          'ind_n_base_part' = 0,
                          'n_base_part' = nbase,
                          'alpha_hi' = 0.1,
                          'beta_hi' = 30,
                          'mu_lat_hi' = 50,
                          'var_lat_lo' = 0.1,
                          'var_lat_hi' = 50,
                          'c_hi' = 100,
                          'd_hi' = 100,
                          'k_1_lo' = 0,
                          'k_1_hi' = 10,
                          'mu_1_hi' = 44,
                          'mu_2_hi' = 0.1,
                          'p_ber_hi' = 1.0,
                          'phi_inf1_hi' = 500,
                          'phi_inf2_hi' = 500,
                          'rho_susc1_hi' = 500,
                          'rho_susc2_hi' = 500,
                          'nu_inf_lo' = 0,
                          'nu_inf_hi' = 1,
                          'tau_susc_lo' = 0,
                          'tau_susc_hi' = 1,
                          'beta_m_hi' = 5,
                          'trace_window' = 20)


## ----inf_setup6--------------------------------------------------------------------------------
para.sf <- data.frame('alpha_sf' = para.key$alpha*1.5,
                      'beta_sf' = para.key$beta*0.55,
                      'mu_lat_sf'	= 2,
                      'var_lat_sf' = 1,
                      'c_sf' = 1.25,
                      'd_sf' = 0.75,
                      'k_1_sf' = 0.01,
                      'mu_1_sf' = 0.2,
                      'mu_2_sf' = 2.5e-6,
                      'p_ber_sf' = 0.02,
                      'phi_inf1_sf' = 1.75,
                      'phi_inf2_sf' = 1.5,
                      'rho_susc1_sf' = 1,
                      'rho_susc2_sf' = 1.25,
                      'nu_inf_sf' = 0.25,
                      'tau_susc_sf' = 0.25,
                      'beta_m_sf' = 1)


## ----inf_setup7--------------------------------------------------------------------------------
t_sample <- sim.out$t_sample

#which individuals are missing genomes and which were sampled
no.genome <- which(t_sample == pars.aux$unassigned_time)
genome.sampled <- which(t_sample != pars.aux$unassigned_time)

print(getwd())
who_sample_list <- c(1,sample(genome.sampled[2:length(genome.sampled)]))
sample_filename <- paste0("sample_list_",subsampling_pct*100,".csv")

if(subsampling_pct < 1){
  print("Random sample, not subset")
  who_sample <- c(1,sample(genome.sampled,size=round(subsampling_pct*length(genome.sampled))-1))
  
  # print("Subset subsample")
  # who_sample <- who_sample_list[1:(round(subsampling_pct*length(genome.sampled)))]
  
  who_no_sample <- setdiff(genome.sampled,who_sample)
  
  write.csv(who_sample,file=sample_filename)
  
  no.genome <- union(which(t_sample == pars.aux$unassigned_time),who_no_sample)
  genome.sampled <- intersect(which(t_sample != pars.aux$unassigned_time),who_sample)
  
  print("unsampled")
  print(no.genome)
  print("sampled")
  print(genome.sampled)
  print("exact percent sampled")
  print(length(genome.sampled)/length(which(t_sample != pars.aux$unassigned_time)))
} else{
  print(who_sample_list)
  write.csv(who_sample_list,file=sample_filename)
}

#the sampling times of those with genomes
round(t_sample[-no.genome], 3)

#setup a matrix to hold 1 sequence per individual or missing data
seq_mat<-matrix(nrow=nrow(epi.data), ncol=pars.aux$n_base)

#identify and store the data from sequences at the time of sampling
#the next step takes approximately 30 seconds to run
for(i in 1:nrow(epi.data)){
  k <- i-1
  if(i %in% no.genome){
    #these are not used, given t.sample for these = the unassigned time
    #n denotes any base (in the IUPAC system)
    seq_mat[i,]<-rep("n", pars.aux$n_base)
  }else{
    #identify the closest corresponding sample to the time of
    #sampling for each sampled individual
    fn<-paste0(sim.out$outputPath,"subject_",k,"_t_nt.txt")
    ts<-as.numeric(unlist(read.csv(fn, header=F)))
    tsamp <- which.min(abs(ts - t_sample[i]))
    #read in the corresponding nucleotide sequence
    fn<-paste0(sim.out$outputPath,"subject_",k,"_nt.txt")
    nts<-read.csv(fn, header = F, stringsAsFactors = F)
    #store the closest corresponding sample for this individual
    seq_mat[i,]<-as.numeric(nts[tsamp,])
  }
  if((i%%2) == 0) {cat(".",sep=""); flush.console()}  #
}


## ----inf_setup8--------------------------------------------------------------------------------

#convert the sequence data to their nucleotide base letters
seq_mat[seq_mat==1]<-"a"
seq_mat[seq_mat==2]<-"g"
seq_mat[seq_mat==3]<-"t"
seq_mat[seq_mat==4]<-"c"

#write row.names to the sequence matrix object, including the 'k'
#and time of sampling. These will be used as tip labels
row.names(seq_mat)<- paste0(epi.data$k, "_", round(sim.out$t_sample,0))


if(!dir.exists("../Infer_Test/gen_inputs")) {
  dir.create("../Infer_Test/gen_inputs")
}


#write these out to a fasta file in the dnaPath directory
library(ape)
write.dna(seq_mat, file='../Infer_Test/gen_inputs/seqs.fasta', format = 'fasta', nbcol = -1)

# if(ndiff_known){
#   filenames <- list.files(path="./outputs/")
#   
#   nbpdiff <- list()
#   dtdiff <- list()
#   tseq <- list()
#   
#   for(i in 0:(npop-1)){
#     ntname <- paste0("./outputs/subject_",i,"_nt.txt")
#     tntname <- paste0("./outputs/subject_",i,"_t_nt.txt")
#     
#     file_lines <- readLines(ntname)
#     tnt_lines <- readLines(tntname)
#     
#     if(length(file_lines) > 0 & length(tnt_lines) > 0){
#       seqlist <- read.csv(ntname,header=F)
#       nseq <- nrow(seqlist)
#       if(nseq > 1){
#         ndiff <- rowSums(seqlist[1:(nseq-1),] != seqlist[2:(nseq),])
#         tlist <- read.csv(tntname,header=F)
#         tdiff <- tlist[2:nseq,1] - tlist[1:(nseq-1),1]
#         
#         nbpdiff[[i+1]] <- ndiff
#         dtdiff[[i+1]] <- tdiff
#         tseq[[i+1]] <- tlist[,1]
#       } else{
#         nbpdiff[[i+1]] <- c(-1)
#         dtdiff[[i+1]] <- c(-1)
#         tseq[[i+1]] <- c(-1)
#       }
#     } else{
#       nbpdiff[[i+1]] <- c(-1)
#       dtdiff[[i+1]] <- c(-1)
#       tseq[[i+1]] <- c(-1)
#     }
#   }
#   
#   if(!file.exists("../Infer_Test/nbpdiff_true.txt")){
#     file.create("../Infer_Test/nbpdiff_true.txt")
#   } else{
#     unlink("../Infer_Test/nbpdiff_true.txt")
#     file.create("../Infer_Test/nbpdiff_true.txt")
#   }
#   
#   fileconn <- "../Infer_Test/nbpdiff_true.txt"
#   
#   # cat(as.character(0),file=fileconn,sep="")
#   # cat(as.character("\n"),file=fileconn,append=T)
#   cat(as.character(nbpdiff[[1]]),file=fileconn,sep=",",append=T)
#   
#   for(l in 2:npop){
#     
#     if(is.null(nbpdiff[[l]])){
#       next
#     }
#     
#     # if((nbpdiff[[l]][1]) >= 0){
#     # cat(as.character(l-1),file=fileconn,sep="",append=T)
#     cat(as.character("\n"),file=fileconn,append=T)
#     cat(as.character(nbpdiff[[l]]),file=fileconn,sep=",",append=T)
#     # }
#   }
#   
#   if(!file.exists("../Infer_Test/t_true.txt")){
#     file.create("../Infer_Test/t_true.txt")
#   } else{
#     unlink("../Infer_Test/t_true.txt")
#     file.create("../Infer_Test/t_true.txt")
#   }
#   
#   fileconn <- "../Infer_Test/t_true.txt"
#   
#   # cat(as.character(0),file=fileconn,sep="")
#   # cat(as.character("\n"),file=fileconn,append=T)
#   cat(as.character(tseq[[1]]),file=fileconn,sep=",",append=T)
#   
#   for(l in 2:npop){
#     if(is.null(tseq[[l]])){
#       next
#     }
#     
#     # if((tseq[[l]][1]) >= 0){
#     # cat(as.character(l-1),file=fileconn,sep="",append=T)
#     # cat(as.character("\n"),file=fileconn,append=T)
#     cat(as.character("\n"),file=fileconn,append=T)
#     cat(as.character(tseq[[l]]),file=fileconn,sep=",",append=T)
#     # }
#   }
# }

## ----inf_run1, eval=F--------------------------------------------------------------------------
pars.aux$n_iterations = niter
pars.aux$n_frequ = nfrequ
pars.aux$n_output_source = noutputsource
pars.aux$n_output_gm = noutputgm
pars.aux$n_cout = ncout

run_params$expected_lambda <- (para.key$mu_1 + para.key$mu_2)*pars.aux$n_base
run_params$nburn <- nburn

print(para.key.inits$beta)

if(subsampling_pct<1){
  sim.out$t_sample[no.genome] <- 9e6
}

save(run_params,
     pars.aux,
     para.key,
     para.key.inits,
     para.priors,
     para.sf,
     file="../run_settings.RData")

setwd("../Infer_Test/")

sink(file="borisout.txt")
time0 <- Sys.time()
infer.out<-infer(covariates = epi.data,
                 moves.inputs = moves.inputs,
                 parsAux = pars.aux,
                 keyInits = para.key.inits,
                 priors = para.priors,
                 scalingFactors = para.sf,
                 seed = seed+1,
                 accTable = sim.out$infected_source,
                 t.sample = sim.out$t_sample,
                 inputPath = "./inputs",
                 outputPath = "./outputs",
                 dnaPath = "./gen_inputs",
                 dnaReference = F)
time1 <- Sys.time()

print(time1-time0)

sink(file=NULL)

speed_values <- data.frame(run_speed=difftime(time0,time1,units="mins"),format="minutes",npop=npop,niter=pars.aux$n_iterations)

save(speed_values,file="speed_test_results.RData") 

if(!cluster){
  source("~/Documents/Research/Phylodynamics/BORIS2/test/parameter_capture.R")
  sink(file="posterior_tree_capture.txt")
  source("~/Documents/Research/Phylodynamics/BORIS2/test/posterior_tree_capture.R")
  sink(NULL)
  source("~/Documents/Research/Phylodynamics/BORIS2/test/bubble_plot.R")
} else{
  ### Calculate parameter traceplots
  ## set to Infer_Test directory if running individually
  raw_params <- read.table("./outputs/parameters_current.log",header=T)
  true_params <- read.csv("../Sim_Test/inputs/parameters_key.csv",header=T)
  t_sample <- read.csv("./inputs/t_sample.csv",header=F)[,1]
  
  load("../run_settings.RData")
  
  epi_sim <- read.csv("../Sim_Test/outputs/epi_sim.csv")
  n_infect_observed <- min(c(sum(epi_sim$t_e <= pars.aux$t_max),sum(t_sample <= pars.aux$t_max)))
  
  # params_current <- raw_params[run_params$nburn:nrow(raw_params),]
  
  params_current <- raw_params[(nrow(raw_params)*2/3):nrow(raw_params),]
  
  params_current <- raw_params
  
  pdf(paste0("params_posterior_seed",run_params$seed,".pdf"),width=15,height=15)
  par(mfrow=c(3,3))
  
  hist(params_current$alpha,main=paste0("alpha = ",true_params$alpha),xlab="")
  abline(v=true_params$alpha,col="blue")
  
  hist(params_current$beta,main=paste0("beta = ",true_params$beta),xlab="")
  abline(v=true_params$beta,col="blue")
  
  hist(params_current$k_1,main=paste0("kappa = ",true_params$k_1),xlab="")
  abline(v=true_params$k_1,col="blue")
  
  hist(params_current$lat_mu,main=paste0("lat_mu = ",true_params$a),xlab="")
  abline(v=true_params$a,col="blue")
  
  hist(params_current$lat_sd,main=paste0("lat_sd = ",true_params$b),xlab="")
  abline(v=true_params$b,col="blue")
  
  hist(params_current$c,main=paste0("c = ",true_params$c),xlab="")
  abline(v=true_params$c,col="blue")
  
  hist(params_current$d,main=paste0("d = ",true_params$d),xlab="")
  abline(v=true_params$d,col="blue")
  
  if("lambda" %in% names(params_current)){
    
    nbpdiffs <- c()
    deltats <- c()
    
    # if(file.exists("nbpdiff_true.txt")){
    #   
    #   nbpdiff_all <- read.csv("nbpdiff_true.txt",header=F,col.names = 1:20)
    #   infect <- which(nbpdiff_all[,1] >= 0)
    #   nbpdiffs <- nbpdiff_all[infect,]
    #   
    #   t_all <- read.csv("t_true.txt",header=F,col.names = 1:20)
    #   infect <- which(t_all[,1] >= 0)
    #   ts <- t_all[infect,]
    #   deltats <- apply(ts,1,diff)
    #   
    #   lambda_mle <- sum(nbpdiffs,na.rm=T)/sum(deltats,na.rm=T)
    # }
    
    hist(params_current$lambda,main=paste0("lambda"),xlab="")
    abline(v=(true_params$mu_1+true_params$mu_2)*run_params$nbase,col="blue",lty=2)
    # abline(v=lambda_mle,col="blue",lty=1)
    
  } else{
    hist(params_current$mu_1,main="Mu_1",xlab="")
    abline(v=true_params$mu_1)
    
    hist(params_current$mu_2,main="Mu_2",xlab="")
    abline(v=true_params$mu_2)
  }
  
  par(mfrow=c(1,1))
  
  par(mfrow=c(3,3))
  
  plot(params_current$alpha,main=paste0("alpha = ",true_params$alpha),xlab="",type="l")
  abline(h=true_params$alpha,col="blue")
  
  plot(params_current$beta,main=paste0("beta = ",true_params$beta),xlab="",type="l")
  abline(h=true_params$beta,col="blue")
  
  plot(params_current$k_1,main=paste0("kappa = ",true_params$k_1),xlab="",type="l")
  abline(h=true_params$k_1,col="blue")
  
  plot(params_current$lat_mu,main=paste0("lat_mu = ",true_params$a),xlab="",type="l")
  abline(h=true_params$a,col="blue")
  
  plot(params_current$lat_sd,main=paste0("lat_sd = ",true_params$b),xlab="",type="l")
  abline(h=true_params$b,col="blue")
  
  plot(params_current$c,main=paste0("c = ",true_params$c),xlab="",type="l")
  abline(h=true_params$c,col="blue")
  
  plot(params_current$d,main=paste0("d = ",true_params$d),xlab="",type="l")
  abline(h=true_params$d,col="blue")
  
  if("lambda" %in% names(params_current)){
    plot(params_current$lambda,main=paste0("lambda"),xlab="",type="l")
    abline(h=(true_params$mu_1+true_params$mu_2)*run_params$nbase,col="blue",lty=2)
    # abline(h=lambda_mle,col="blue",lty=1)
  } else{
    plot(params_current$mu_1,main=paste0("Mu_1"),xlab="",type="l")
    abline(h=true_params$mu_1,col="blue",lty=2)
    
    plot(params_current$mu_2,main=paste0("Mu_2"),xlab="",type="l")
    abline(h=true_params$mu_2,col="blue",lty=2)
  }
  
  par(mfrow=c(1,1))
  
  plot(params_current$log_likelihood,type="l",xlab="Iteration",main="Log-Likelihood")
  
  dev.off()
  
  t_e_raw <- read.csv("./outputs/t_e_current.csv",header=F)
  
  run_params$nburn = 0
  
  t_e_current <- t_e_raw[(nrow(t_e_raw)-(nrow(t_e_raw)-run_params$nburn)):nrow(t_e_raw),]
  
  pdf(file=paste0("exposure_posterior_seed",run_params$seed,".pdf"),width = 20,height=20)
  
  par(mfrow=c(5,5))
  
  for(i in 1:ncol(t_e_current)){
    hist(t_e_current[,i],main=i-1)
    abline(v=epi_sim$t_e[i],col="blue")
  }
  
  dev.off()
  
  pdf(file=paste0("exposure_trace_seed",run_params$seed,".pdf"),width = 20,height=20)
  
  par(mfrow=c(5,5))
  
  for(i in 1:ncol(t_e_current)){
    plot(t_e_current[,i],main=i-1,type="l")
    abline(h=epi_sim$t_e[i],col="blue")
  }
  
  dev.off()
  
  png(filename="params.png",width = 900,height = 400)
  par(mfrow=c(2,4))
  hist(params_current$alpha,main=expression(alpha),xlab="")
  abline(v=true_params$alpha,col="blue")
  plot(params_current$alpha,xlab="",type="l",main=expression(alpha),ylab="")
  abline(h=true_params$alpha,col="blue")
  
  hist(params_current$beta,main=expression(beta),xlab="")
  abline(v=true_params$beta,col="blue")
  plot(params_current$beta,xlab="",type="l",main=expression(beta),ylab="")
  abline(h=true_params$beta,col="blue")
  
  hist(params_current$k_1,main=expression(kappa),xlab="")
  abline(v=true_params$k_1,col="blue")
  plot(params_current$k_1,xlab="",type="l",main=expression(kappa),ylab="")
  abline(h=true_params$k_1,col="blue")
  
  if("lambda" %in% names(params_current)){
    hist(params_current$lambda,main=expression(lambda),xlab="")
    abline(v=(true_params$mu_1+true_params$mu_2)*run_params$nbase,col="blue",lty=2)
    plot(params_current$lambda,xlab="",type="l",main=expression(lambda),ylab="")
    abline(h=(true_params$mu_1+true_params$mu_2)*run_params$nbase,col="blue",lty=2)
  } else{
    
  }
  par(mfrow=c(1,1))
  dev.off()
  
}
