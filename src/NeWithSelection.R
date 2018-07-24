
source("src/fun/ColorBlindPalette.R")



load(file="results/simtable.RData")



pdf(file="results/NeWithSelection.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(49,50,51,52,1)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  sel_coef <- sim_table$sel_coef[scenarios[i]]
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))

  
  if (i==1){
    
    boxplot( results$Ne_hat/2,
             at = sel_coef,
             outline = FALSE,
             boxwex = 0.05,
             notch = T,
             col = 1,
             border = 1,
             names = c("0"),
             ylab = expression("effective population size, "*italic(N)[e]),
             xlab = expression("selection coefficient, "*italic(s)),
             cex.lab=0.75,
             ylim = c(0,1000),
             xlim = c(0.1,0.5) ,cex.axis=0.7,las=1)
    
  }else{
    boxplot( results$Ne_hat/2,
             at = sel_coef,
             outline = FALSE,
             axes = FALSE,
             boxwex = 0.05,
             notch = T,
             col = 1,
             border = 1,
             add=T)
    
  }
  
}
axis(side=1,at=c(0.1,0.2,0.3,0.4,0.5),las=2,cex.axis=0.6)
abline(h=500,lty=2)

legend(x=1,y=1,
       legend=c(expression("s="* 0.1),
                expression("s="* 0.2),
                expression("s="* 0.3),
                expression("s="* 0.4),
                expression("s="* 0.5)),
       lty=1,
       col=1:5,
       cex=0.75,
       bty="n")

dev.off()

