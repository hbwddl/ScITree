## ----setup--------------------------------------------------------------------
## Necessary packages
library(ScTree)

library(ggplot2)
library(ape)
library(dplyr)

base_dir <- "~/Documents/Research/Phylodynamics/ScITree/vignettes/"

setwd(base_dir)

## -----------------------------------------------------------------------------
## Run setup
## Simulation seed
random_seed <- 108

set.seed(seed=random_seed)

## -----------------------------------------------------------------------------
n_pop <- 150 # Number of individuals in the population
nbase <- 8000 # Number of pathogen nucleotide bases
tmax <- 100 # Maximum time of outbreak

alphaval <- 0.0001 # Background infection rate
betaval <- 8 # Infectious pressure
kappaval <- 0.02 # Exponential kernel decay rate
aval <- 10 # Exposed sojourn shape (gamma distribution)
bval <- 0.5 # Exposed sojourn scale (gamma distribution)
cval <- 2 # Infectious sojourn shape (Weibull distribution)
dval <- 2 # Infectious sojourn scale (Weibull distribution)
mu1val <- 0.002 # Kimura substitution model transition rate
mu2val <- 0.0005 # Kimura substitution model transversion rate
pval <- 0.0001 # Mutation probability from "grand master" background sequence

maxx <- 2000 # Maximum x coordinate
maxy <- 2000 # Maximum y coordinate

## -----------------------------------------------------------------------------
sim.epi.data <- data.frame(k=0:(n_pop-1), ## Index
                           coor_x=runif(n_pop,0,maxx), ## X coordinates
                           coor_y=runif(n_pop,0,maxy), ## Y coordinates
                           t_e=c(0,rep(9e+06,n_pop-1)), ## Exposure time
                           t_i=rep(9e+06,n_pop), ## Infectious time
                           t_r=rep(9e+06,n_pop), ## Removal time
                           ftype0=sample(c(0,1),n_pop,replace=T), ## Unused
                           herdn=sample(min(sim.epi.input$herdn):max(sim.epi.input$herdn), 
                                        n_pop,
                                        replace=T), ## Unused
                           status=c(2,rep(1,n_pop-1))) ## Unused

sim.epi.data$ftype1=ifelse(sim.epi.data$ftype0==1,
              0,
              sample(c(0,1),n_pop,replace=T)) ## Unused
sim.epi.data$ftype2=ifelse(sim.epi.data$ftype0==0 & sim.epi.data$ftype1==0,
              1,
              0) ## Unused

data(sim.epi.input)
sim.epi.data <- sim.epi.data[,names(sim.epi.input)]

## -----------------------------------------------------------------------------
niter = 20 # Number of iterations
nfrequ = round(niter/2) # Number of individuals to update at each iteration
noutputsource = max(round(niter/10000),1) # Frequency of source output
noutputgm = max(round(niter/1000),1) # Frequency of grand master sequence output (not used)
ncout = max(round(niter/1000),1) # Frequency of console output

nburn = round(niter/10) # Number of burnin iterations

## -----------------------------------------------------------------------------
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

pars.aux <- data.frame('n' = nrow(sim.epi.data),
                       'seed' = random_seed,
                       'n_base' = nbase,
                       'n_seq' = 5,
                       't_max' = tmax,
                       'unassigned_time' =  9e+6,
                       'sample_range' = 10,
                       'partial_seq_out' = 0,
                       'n_base_part' = nbase,
                       'n_index' = 1,
                       'coord_type' = 'cartesian',
                       'kernel_type' = 'exponential',
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

## -----------------------------------------------------------------------------
sim.out <- sim(epi.inputs = sim.epi.data,
               moves.inputs = sim.moves.input,
               parsKey = para.key,
               parsAux = pars.aux,
               inputPath = "./inputs",
               outputPath = "./outputs")


## -----------------------------------------------------------------------------
arrow_outbreak_plot <- function(epi.data,true.source.c,t.max=100){
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
           guides(colour="none",is.infected="none",alpha="none") +
           xlab("X") + ylab("Y"))
}

outbreak_plot <- arrow_outbreak_plot(sim.out$epi.sim,sim.out$infected_source,tmax)
print(outbreak_plot)

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


## -----------------------------------------------------------------------------
dir.create("../Infer_Test")

fn <- paste0(sim.out$outputPath, "subject_0_nt.txt")
seqs1 <- as.matrix(read.csv(fn, header=F))

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

write.dna(seqs1, file='seqs1.fasta', format = 'fasta', nbcol = -1)

epi.data <- as.data.frame(sim.out$epi.sim)

#make sure indexing is correct
epi.data$k<-1:nrow(epi.data)-1

#produce the required initial values for the temporal fields.
#assume that exposure has occurred 24 hours before the onset of infectiousness:
epi.data$t_e <- epi.data$t_i - 1

epi.data$t_e[epi.data$t_e > 10000] <- pars.aux$unassigned_time
epi.data$t_r <- round(epi.data$t_r, 0)

# UNUSED
#wrangle the farm type variables into the required format
epi.data$ftype <- NA
epi.data$ftype[epi.data$ftype0 == 1] <- 0
epi.data$ftype[epi.data$ftype1 == 1] <- 1
epi.data$ftype[epi.data$ftype2 == 1] <- 2

epi.data<-epi.data[,c('k','coor_x','coor_y','t_e','t_i','t_r','ftype','herdn')]
head(epi.data)

## -----------------------------------------------------------------------------
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

## UNUSED
#use the movement data from previously
moves.inputs<-sim.out$moves.inputs
head(moves.inputs)

## -----------------------------------------------------------------------------
pars.aux <- data.frame('n' = nrow(epi.data),
                       'kernel_type' = 'exponential',
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

para.key.inits <- data.frame('alpha' = para.key$alpha,
                             'beta' = para.key$beta/2,
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

## -----------------------------------------------------------------------------
t_sample <- sim.out$t_sample

no.genome <- which(t_sample == pars.aux$unassigned_time)
genome.sampled <- which(t_sample != pars.aux$unassigned_time)

#identify and store the data from sequences at the time of sampling
#the next step takes approximately 30 seconds to run
seq_mat<-matrix(nrow=nrow(epi.data), ncol=pars.aux$n_base)

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

write.dna(seq_mat, file='../Infer_Test/gen_inputs/seqs.fasta', format = 'fasta', nbcol = -1)



## -----------------------------------------------------------------------------
setwd(base_dir)
setwd("./Infer_Test/")

## ----results='hide'-----------------------------------------------------------
pars.aux$n_iterations = niter
pars.aux$n_frequ = nfrequ
pars.aux$n_output_source = noutputsource
pars.aux$n_output_gm = noutputgm
pars.aux$n_cout = ncout

sink(file = "inference_output.txt")
infer.out<-infer(covariates = epi.data,
                 moves.inputs = moves.inputs,
                 parsAux = pars.aux,
                 keyInits = para.key.inits,
                 priors = para.priors,
                 scalingFactors = para.sf,
                 seed = random_seed + 1,
                 accTable = sim.out$infected_source,
                 t.sample = sim.out$t_sample,
                 inputPath = "./inputs",
                 outputPath = "./outputs",
                 dnaPath = "./gen_inputs",
                 dnaReference = F)
sink(file=NULL)

## -----------------------------------------------------------------------------
source("../../test/parameter_capture.R")
source("../../test/posterior_tree_capture.R")
source("../../test/bubble_plot.R")

