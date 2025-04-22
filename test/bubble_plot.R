## Bubble plot
## Set to Infer_test directory
library(ggplot2)
library(dplyr)

inf_source <- read.csv("./outputs/infected_source_current.csv",header=F)

if(dir.exists("../Sim_Test/outputs/")){
  true_source <- as.numeric(read.csv("../Sim_Test/outputs/infected_source.txt",header=F)[,1])
} else{
  true_source <- as.numeric(read.csv("./inputs/atab_from.csv",header=F)[,1])
}

epi_data = read.csv("./inputs/epi.csv",header=T)

true_source_9999 <- true_source
true_source_9999[true_source == -99] <- 9999

niter <- nrow(inf_source)

inf_source_burnin <- inf_source[round(niter/10):nrow(inf_source),]

inf_source_burnin <- inf_source

pct_correct_source <- true_source
max_post_prob <- true_source
second_post_prob <- true_source
pct_top2 <- cbind(true_source,true_source)

for(i in 1:length(max_post_prob)){
  post_dist <- table(inf_source_burnin[,i])
  # print(post_dist)
  max_post_prob[i] <- as.numeric(names(post_dist)[which.max(post_dist)])
  sort_posterior_src <- order(post_dist,decreasing=T)

  sort_post_dist <- sort(post_dist,decreasing=T)

  pct_top2[i,] <- sort_post_dist[c(1,2)]/sum(sort_post_dist)

  second_post_prob[i] <- as.numeric(names(post_dist)[sort_posterior_src[2]])

  which_correct_source <- which(as.numeric(names(post_dist)) == true_source_9999[i])

  if(length(which_correct_source) > 0){
    pct_correct_source[i] = post_dist[which_correct_source]/sum(post_dist);
  } else{
    pct_correct_source[i] = 0
  }
}

epi_data$pct_correct = pct_correct_source
epi_data$is_infected = ifelse(true_source >= 0,"yes","no")
epi_data$true_source = true_source

## Make the arrow start and ends for everyone
true_source_r <- true_source + 1
true_source_r[true_source_r > nrow(epi_data) | true_source_r < 0] <- NA
epi_data$arrow_start_x <- epi_data$coor_x[true_source_r]
epi_data$arrow_start_y <- epi_data$coor_y[true_source_r]

bubble_plot <- epi_data %>% ggplot(aes(x=coor_x,y=coor_y,alpha=0.83)) +
  geom_point(mapping=aes(colour=is_infected)) +
  scale_color_manual(values=c("no"="grey85","yes"="white")) +
  geom_point(data=subset(epi_data,is_infected=="yes"),mapping=aes(size=pct_correct)) +
  scale_size_continuous(range=c(0.5,3)) +
  theme_classic() +
  geom_point(data=subset(epi_data,true_source==9999),shape=3) +
  labs(size="Percent Correct") +
  guides(colour="none",is_infected="none",alpha="none") +
  xlab("X") +
  ylab("Y")

bubble_plot_arrow <- epi_data %>% ggplot(aes(x=coor_x,y=coor_y,alpha=0.83)) +
  geom_point(mapping=aes(colour=is_infected)) +
  scale_color_manual(values=c("no"="grey85","yes"="white")) +
  geom_point(data=subset(epi_data,is_infected=="yes")) +
  # geom_point(data=subset(epi_data,is_infected=="yes"),mapping=aes(size=pct_correct)) +
  # scale_size_continuous(range=c(0.5,3)) +
  theme_classic() +
  geom_point(data=subset(epi_data,true_source==9999),shape=3) +
  # labs(size="Percent Correct") +
  labs(linewidth="Percent Correct") +
  geom_segment(aes(x=arrow_start_x,y=arrow_start_y,xend=coor_x,yend=coor_y,linewidth=pct_correct),arrow=arrow(length=unit(0.2,'cm')),color="darkblue") +
  scale_linewidth_continuous(range=c(0.25,2)) +
  guides(colour="none",is_infected="none",alpha="none") +
  xlab("X") +
  ylab("Y")

pdf(file="./bubble_plot.pdf",width=8,height=8)
print(bubble_plot)
print(bubble_plot_arrow)
dev.off()


bubble_plot <- epi_data %>% ggplot(aes(x=coor_x,y=coor_y,alpha=0.83)) +
  xlim(144.825,145) +
  ylim(-36.475,-36.35) +
  geom_point(mapping=aes(colour=is_infected)) +
  scale_color_manual(values=c("no"="grey85","yes"="white")) +
  geom_point(data=subset(epi_data,is_infected=="yes"),mapping=aes(size=pct_correct)) +
  scale_size_continuous(range=c(0.5,3)) +
  theme_classic() +
  geom_point(data=subset(epi_data,true_source==9999),shape=3) +
  labs(size="Percent Correct") +
  guides(colour="none",is_infected="none",alpha="none") +
  xlab("X") +
  ylab("Y")

bubble_plot_arrow <- epi_data %>% ggplot(aes(x=coor_x,y=coor_y,alpha=0.83)) +
  geom_point(mapping=aes(colour=is_infected)) +
  scale_color_manual(values=c("no"="grey85","yes"="white")) +
  geom_point(data=subset(epi_data,is_infected=="yes")) +
  # geom_point(data=subset(epi_data,is_infected=="yes"),mapping=aes(size=pct_correct)) +
  # scale_size_continuous(range=c(0.5,3)) +
  theme_classic() +
  geom_point(data=subset(epi_data,true_source==9999),shape=3) +
  # labs(size="Percent Correct") +
  labs(linewidth="Percent Correct") +
  geom_segment(aes(x=arrow_start_x,y=arrow_start_y,xend=coor_x,yend=coor_y,linewidth=pct_correct),arrow=arrow(length=unit(0.15,'cm')),color="darkblue") +
  scale_linewidth_continuous(range=c(0.25,2)) +
  guides(colour="none",is_infected="none",alpha="none") +
  xlab("X") +
  ylab("Y")

hist(pct_correct)


