### Calculate parameter traceplots
## set to Infer_Test directory if running individually
library(coda)

raw_params <- read.table("./outputs/parameters_current.log",header=T)
true_params <- read.csv("../Sim_Test/inputs/parameters_key.csv",header=T)
aux_params <- read.csv("../Sim_Test/inputs/parameters_other.csv",header=T)
t_sample <- read.csv("./inputs/t_sample.csv",header=F)[,1]

epi_sim <- read.csv("../Sim_Test/outputs/epi_sim.csv")
n_infect_observed <- min(c(sum(epi_sim$t_e <= pars.aux$t_max),sum(t_sample <= pars.aux$t_max)))

params_current <- raw_params[(nrow(raw_params)*2/3):nrow(raw_params),]

params_current <- raw_params

params.mcmc <- as.mcmc(true_params)

pdf("params_posterior_plots.pdf",width=15,height=15)
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
  
    if(file.exists("nbpdiff_true.txt")){
      
      nbpdiff_all <- read.csv("nbpdiff_true.txt",header=F,col.names = 1:20)
      infect <- which(nbpdiff_all[,1] >= 0)
      nbpdiffs <- nbpdiff_all[infect,]
      
      t_all <- read.csv("t_true.txt",header=F,col.names = 1:20)
      infect <- which(t_all[,1] >= 0)
      ts <- t_all[infect,]
      deltats <- apply(ts,1,diff)
      
      lambda_mle <- sum(nbpdiffs,na.rm=T)/sum(deltats,na.rm=T)
    } else{
      lambda_mle <- (true_params$mu_1+2*true_params$mu_2)*aux_params$n_base
    }
    
    hist(params_current$lambda,main=paste0("lambda"),xlab="")
    abline(v=(true_params$mu_1+2*true_params$mu_2)*aux_params$n_base,col="blue",lty=2)
    abline(v=lambda_mle,col="blue",lty=1)

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
  abline(h=(true_params$mu_1+2*true_params$mu_2)*aux_params$n_base,col="blue",lty=2)
  abline(h=lambda_mle,col="blue",lty=1)
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

n_burnin = 0

t_e_current <- t_e_raw[(nrow(t_e_raw)-(nrow(t_e_raw)-n_burnin)):nrow(t_e_raw),]

pdf(file="exposure_posterior_seed.pdf",width = 20,height=20)

par(mfrow=c(5,5))

for(i in 1:ncol(t_e_current)){
  hist(t_e_current[,i],main=i-1)
  abline(v=epi_sim$t_e[i],col="blue")
}

dev.off()

pdf(file="exposure_traceplots.pdf",width = 20,height=20)

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
  abline(v=(true_params$mu_1+2*true_params$mu_2)*aux_params$n_base,col="blue",lty=2)
  plot(params_current$lambda,xlab="",type="l",main=expression(lambda),ylab="")
  abline(h=(true_params$mu_1+2*true_params$mu_2)*aux_params$n_base,col="blue",lty=2)
} else{
  
}
par(mfrow=c(1,1))
dev.off()

### Acceptance Probability
p_accept_lambda <- sum(params_current$lambda[1:(length(params_current$lambda)-1)] != params_current$lambda[2:(length(params_current$lambda))])/(length(params_current$lambda)-1)
print(paste0("Accept P Lambda: ", p_accept_lambda))
