source("src/fun/ColorBlindPalette.R")
library(dplyr)
library(ggplot2)
library(gridExtra)
load(file="results/simtable.RData")


sigma_values <- c(0.500,0.750,0.800,0.850,
                  0.900,0.950,0.975,1.000)
scenarios <- c(22:26,28:30)
number_of_replicates=100
temp_res1 <- as.data.frame(matrix(NA,nrow=number_of_replicates,ncol=length(scenarios)))
for (i in seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  sigma <- sim_table$sigma[scenarios[i]]
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  temp_res1[,i]<-results$NeHatFst/2
  colnames(temp_res1)[i]<-sigma
}  



scenarios <- c(21,49,50,51,52,1)

number_of_replicates=100
temp_res2 <- as.data.frame(matrix(NA,nrow=number_of_replicates,ncol=length(scenarios)))
for (i in seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  sel_coef <- sim_table$sel_coef[scenarios[i]]
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))

  if (scenarios[i]==21){
    temp_res2[,i]<-results$NeHatFst/2
    colnames(temp_res2)[i]<-sel_coef
  }else{
    temp_res2[,i]<-results$Ne_hat/2
    colnames(temp_res2)[i]<-sel_coef
  }
  
  
  
}


pdf(file="results/Ne.pdf",width=4,height=6)
layout(matrix(c(1,2), 2,1,byrow = TRUE))
par(mar=c(4,4,0.2,0.2)+0.1)

boxplot( temp_res1,
         at = sigma_values,
         outline = FALSE,
         boxwex = 0.02,
         notch = T,
         col = cbPalette1[6],
         border = cbPalette1[6],
         names = NA,
         ylab = expression("effective population size, "*italic(N)[e]),
         xlab = expression("selfing rate, "*sigma),
         cex.lab=0.75,
         ylim = c(0,1000),
         xlim = c(0.5,1)  ,cex.axis=0.7,las=1)


trueNe <- ((2-sigma_values)*500)/2

sigma_values_axis <- sigma_values
#sigma_values_axis[c(7,9)]<-NA

lines(sigma_values,trueNe,lty=2,lwd=1.5)
axis(side=1,at=sigma_values_axis,las=2,cex.axis=0.6)

text(x=0.98,y=950,label="a",cex=1.5)

sel_values <- c(0.0,0.1,0.2,0.3,0.4,0.5)

boxplot( temp_res2,
         at=sel_values,
         outline = FALSE,
         boxwex = 0.03,
         notch = T,
         col = cbPalette1[6],
         border = cbPalette1[6],
         names = NA,
         ylab = expression("effective population size, "*italic(N)[e]),
         xlab = expression("selection coefficient, "*italic(s)),
         cex.lab=0.75,
         ylim = c(0,1000),
         xlim=c(0,0.5),
         cex.axis=0.7,las=1)

N <- 500

sel_values <- c(0.0,0.1,0.2,0.3,0.4,0.5)
#sigma_values_axis[c(7,9)]<-NA

abline(h=N,lty=2,lwd=1.5)
axis(side=1,at=sel_values,las=2,cex.axis=0.6)

text(x=0.49,y=950,label="b",cex=1.5)



dev.off()

#box()



#####################




























data1 <- data.frame(
  sigma=c( rep(0.500,number_of_replicates), 
           rep(0.750,number_of_replicates), 
           rep(0.800,number_of_replicates), 
           rep(0.850,number_of_replicates), 
           rep(0.900,number_of_replicates), 
           rep(0.950,number_of_replicates), 
           rep(0.975,number_of_replicates), 
           rep(1.000,number_of_replicates)  ),
  Ne=c( temp_res1[,1],
        temp_res1[,2],
        temp_res1[,3],
        temp_res1[,4],
        temp_res1[,5],
        temp_res1[,6],
        temp_res1[,7],
        temp_res1[,8] )
)

data2 <- data.frame(
  selection_coef=c( rep(0.0,number_of_replicates), 
                    rep(0.1,number_of_replicates), 
                    rep(0.2,number_of_replicates), 
                    rep(0.3,number_of_replicates), 
                    rep(0.4,number_of_replicates), 
                    rep(0.5,number_of_replicates)  ),
  Ne=c( temp_res2[,1],
           temp_res2[,2],
           temp_res2[,3],
           temp_res2[,4],
           temp_res2[,5],
           temp_res2[,6] )
)

sample_size1 = data1 %>% group_by(sigma) %>% summarize(num=n())
sample_size2 = data2 %>% group_by(selection_coef) %>% summarize(num=n())


trueNe <- ((2-sigma_values)*500)/2


pdf(file="results/NeWithSelection.pdf",width=5,height=8)

p1 <- data1 %>%
  left_join(sample_size1) %>%
  mutate(myaxis = paste0(sigma)) %>%
  ggplot( aes(x=myaxis, y=Ne)) +
  geom_violin(width=1,trim = FALSE) +
  geom_boxplot(width=0.1, color="darkgrey", alpha=0.2) +
  geom_segment(aes(x=0,y=trueNe[1],xend=1.5,yend=trueNe[1]),linetype=2) +
  geom_segment(aes(x=1.5,y=trueNe[2],xend=2.5,yend=trueNe[2]),linetype=2) +
  geom_segment(aes(x=2.5,y=trueNe[3],xend=3.5,yend=trueNe[3]),linetype=2) +
  geom_segment(aes(x=3.5,y=trueNe[4],xend=4.5,yend=trueNe[4]),linetype=2) +
  geom_segment(aes(x=4.5,y=trueNe[5],xend=5.5,yend=trueNe[5]),linetype=2) +
  geom_segment(aes(x=5.5,y=trueNe[6],xend=6.5,yend=trueNe[6]),linetype=2) +
  geom_segment(aes(x=6.5,y=trueNe[7],xend=7.5,yend=trueNe[7]),linetype=2) +
  geom_segment(aes(x=7.5,y=trueNe[8],xend=9,yend=trueNe[8]),linetype=2) +
  #geom_text(aes(label="a", x=8.5,y=1000))+
  annotate("text", x=8.5, y=1000, label= "a",size=8)+
  #geom_hline(yintercept = 500, linetype=2) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab(expression("effective population size, "*italic(N)[e])) +
  xlab(expression("selfing rate, "*sigma)) + 
  ylim(0,1000) +
  #geom_jitter(height = 0, width = 0.2) +
  theme_bw()+ theme_classic() + theme(panel.border = element_rect(linetype = "solid", fill = NA))



p2 <- data2 %>%
  left_join(sample_size2) %>%
  mutate(myaxis = paste0(selection_coef)) %>%
  ggplot( aes(x=myaxis, y=Ne)) +
  geom_violin(width=1,trim = FALSE) +
  geom_boxplot(width=0.1, color="darkgrey", alpha=0.2) +
  geom_hline(yintercept = 500, linetype=2) +
  annotate("text", x=6.3, y=1000, label= "b",size=8)+
  #geom_text(aes(label="b", x=6.3,y=1000))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab(expression("effective population size, "*italic(N)[e])) +
  xlab(expression("selection coefficient, "*italic(s))) + 
  ylim(0,1000) +
  #geom_jitter(height = 0, width = 0.2) +
  theme_bw()+ theme_classic() + theme(panel.border = element_rect(linetype = "solid", fill = NA))

grid.arrange(p1, p2, nrow = 2)
dev.off()

