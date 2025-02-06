### Posterior source summary
library(ggplot2)
library(magrittr)

params_other <- read.csv("./inputs/parameters_other.csv")

setwd("./outputs/")
inf_source <- read.csv("infected_source_current.csv",header=F)

true_source <- as.numeric(read.csv("../inputs/atab_from.csv",header=F)[,1])

true_source_9999 <- true_source
# true_source_9999[true_source == -99] <- 9999

niter <- nrow(inf_source)

inf_source_burnin <- inf_source[round(niter/10):nrow(inf_source),]

inf_source_burnin <- inf_source

max_post_prob <- true_source
second_post_prob <- true_source
pct_top2 <- cbind(true_source,true_source)

pdf(file="../source_posterior.pdf",width=8,height=6,onefile=T)

for(i in 1:length(max_post_prob)){
  post_dist <- table(inf_source_burnin[,i])
  # print(post_dist)
  max_post_prob[i] <- as.numeric(names(post_dist)[which.max(post_dist)])
  sort_posterior_src <- order(post_dist,decreasing=T)

  sort_post_dist <- sort(post_dist,decreasing=T)

  pct_top2[i,] <- sort_post_dist[c(1,2)]/sum(sort_post_dist)

  second_post_prob[i] <- as.numeric(names(post_dist)[sort_posterior_src[2]])

  top5 <- sort_post_dist[1:5]
  
  top5_label <- names(top5)
  
  top5_pct <- as.numeric(top5/sum(sort_post_dist)*100)
  
  is_true <- ifelse(top5_label == true_source[i],"Yes","No")
  
  barplot_dat <- data.frame(top5_label,top5_pct,is_true)
  
  posterior_src_5 <- ggplot(data=barplot_dat,aes(x=factor(top5_label,levels=top5_label),y=top5_pct,fill=is_true)) +
                      geom_bar(stat="identity") +
                      scale_fill_manual(values=c("#999999", "darkgreen")) +
                      labs(fill="Correct Source",x="Source",y="Posterior Probability (%)",title = paste0("Posterior Probability of Source, Individual ",i)) +
                                          theme_minimal()
  
  print(posterior_src_5)
}

dev.off()

all_correct <- sum(max_post_prob == true_source_9999)/length(max_post_prob)

which_infected <- which(true_source != -99)

inf_correct <- sum((max_post_prob == true_source_9999)[which_infected])/length(which_infected)

which_correct_inf <- intersect(which(max_post_prob == true_source_9999),which_infected)

which_incorrect <- which(max_post_prob != true_source_9999)
print("which incorrect (c index)")
print(paste0(which_incorrect-1,collapse=", "))
print("true source")
print(true_source_9999[which_incorrect])
print("estimated source")
print(max_post_prob[which_incorrect])

print("Percent of posterior, first most likely, incorrect")
print(pct_top2[which_incorrect,1])

print("Percent of posterior, second most likely")
print(pct_top2[which_incorrect,2])

print("% correct, highest posterior, all")
print(inf_correct)
print(paste0("Fraction: ",sum((max_post_prob == true_source_9999)[which_infected]),"/",length(which_infected)))

print("% of posterior probabilities, highest posterior, all")
print(pct_top2[which_infected,1])

print("% correct, highest + second highest posterior, all")
two_correct = sum(((max_post_prob == true_source_9999) | (second_post_prob == true_source_9999))[which_infected],na.rm=T)/length(which_infected)
print(two_correct)
print(paste0("Fraction: ",sum(((max_post_prob == true_source_9999) | (second_post_prob == true_source_9999))[which_infected],na.rm=T),"/",length(which_infected)))

print("% of posterior probabilities, second highest posterior, all")
print(pct_top2[which_infected,2])

print("% correct, primary infection")
which_primary = which(true_source == 9999)
print(sum((max_post_prob == true_source_9999)[which_primary])/length(which_primary))
print(paste0("Fraction: ",sum((max_post_prob == true_source_9999)[which_primary]),"/",length(which_primary)))

print("% correct, secondary infection")
which_2ndary = which(true_source >= 0 & true_source != 9999)
print(sum((max_post_prob == true_source_9999)[which_2ndary])/length(which_2ndary))
print(paste0("Fraction: ",sum((max_post_prob == true_source_9999)[which_2ndary]),"/",length(which_2ndary)))

print("% posterior, correct max posterior probability")
print(pct_top2[which_correct_inf,])

print("% posterior, incorrect first correct second")
which_2corr = which((max_post_prob != true_source_9999) & (second_post_prob == true_source_9999))
print(pct_top2[which_2corr,])

# #### Count correct % by iteration, get posterior
is_perfect <- rep(NA,nrow(inf_source_burnin))

n_correct <- rep(NA,nrow(inf_source_burnin))
for(i in 1:nrow(inf_source_burnin)){
  n_correct[i] <- sum(inf_source_burnin[i,which_infected] == true_source[which_infected])
  is_perfect[i] <- sum(inf_source_burnin[i,which_infected] == true_source[which_infected]) == length(which_infected)
}
 
pct_correct <- n_correct/length(which_infected)

pdf(file="../pct_source_correct.pdf",width=10,height=8)

hist(pct_correct,main=paste0("Percent correct each iteration\n% correct identified by posterior probability",inf_correct))
plot(pct_correct,main="Percent correct",xlab="Iteration",ylab="Percent",type="l")

hist(pct_correct,main="Percent correct each iteration")
plot(pct_correct,main="Percent correct",xlab="Iteration",ylab="Percent",type="l")
dev.off()

print("perfect trees")
print(sum(is_perfect)/nrow(inf_source_burnin))

### Source trace
inf_source_negative <- inf_source_burnin
true_source_negative <- true_source_9999

inf_source_negative[inf_source_negative==9999] <- -20
inf_source_negative[true_source_negative==-99] <- -30

true_source_negative[true_source_negative==9999] <- -20
true_source_negative[true_source_negative==-99] <- -30

n_cluster_mcmc <- rowSums(inf_source_negative == -20)

pdf(file="../source_trace.pdf",width=20,height=20)
par(mfrow=c(4,4))

for(i in 1:ncol(inf_source_negative)){
  plot(inf_source_negative[,i],type="l",main=i-1,xlab="Iteration",ylab="Source",ylim=c(-21,ncol(inf_source_negative)+1))
  abline(h=true_source_negative[i],col="blue")
}

par(mfrow=c(1,1))
dev.off()

setwd("..")

epi.data <- read.csv("inputs/epi.csv")

epi.data$premise <- 0:(nrow(epi.data)-1)

source_current <- read.csv("./outputs/infected_source_current.csv",header=F)

source_posterior <- apply(source_current,2,table)
epi.data$source_estimate <- unlist(lapply(source_posterior,function(x) return(as.numeric(names(x)[which.max(x)]))))
epi.data$source_estimate_r <- epi.data$source_estimate+1
epi.data$source_estimate_p <- unlist(lapply(source_posterior,function(x) return(as.numeric(x[which.max(x)]))))/nrow(source_current)

epi.data$source_estimate_r[epi.data$source_estimate_r > nrow(epi.data) | epi.data$source_estimate_r < 0] <- NA
epi.data$arrow_start_x <- epi.data$coor_x[epi.data$source_estimate_r]
epi.data$arrow_start_y <- epi.data$coor_y[epi.data$source_estimate_r]

source_true <- read.csv("inputs/atab_from.csv",header=F)[,1]

epi.data$source_true <- source_true
epi.data$source_true_r <- source_true+1

epi.data$source_true_r[epi.data$source_true_r > nrow(epi.data) | epi.data$source_true_r < 0] <- NA
epi.data$arrow_start_x_true <- epi.data$coor_x[epi.data$source_true_r]
epi.data$arrow_start_y_true <- epi.data$coor_y[epi.data$source_true_r]

n_infected <- sum(source_true >= 0)
n_cluster <- sum(source_true == 9999)

n_cluster_true <- sum(epi.data$source_true == 9999 & epi.data$t_e <= params_other$t_max)
n_cluster_estimate <- sum(epi.data$source_estimate == 9999 & epi.data$t_e <= params_other$t_max)

print("Number of clusters estimates, true/estimated")
print(c(n_cluster_true,n_cluster_estimate))

estimate_plot_arrow <- epi.data %>% ggplot(aes(x=coor_x,y=coor_y)) +
  geom_point() +
  geom_text(aes(label=premise),hjust=0, vjust=0) +
  theme_classic() +
  labs(linewidth="Percent of Posterior") +
  geom_segment(aes(x=arrow_start_x_true,y=arrow_start_y_true,xend=coor_x,yend=coor_y,linewidth=.3),arrow=arrow(length=unit(0.2,'cm')),color="darkred") +
  geom_segment(aes(x=arrow_start_x,y=arrow_start_y,xend=coor_x,yend=coor_y,linewidth=source_estimate_p,alpha=0.7),arrow=arrow(length=unit(0.2,'cm')),color="blue") +
  scale_linewidth_continuous(range=c(0.25,2)) +
  guides(colour="none",is_infected="none",alpha="none")

print(estimate_plot_arrow)

pdf(file="most_probable_tree_accuracy.pdf",width=12,height=8)

print(estimate_plot_arrow)

dev.off()

### Number of clusters
print(paste0("True: ",sum(true_source_9999 == 9999)))
print(paste0("Estimated: ",sum(max_post_prob[which_infected] == 9999)))

png("n_cluster_mcmc.png",width=700,height=700,res=100)
plot(n_cluster_mcmc,type="l",main="Number of clusters each iteration")
dev.off()