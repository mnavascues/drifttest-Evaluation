
source("src/fun/ColorBlindPalette.R")
palette(cbPalette2)

load(file="results/simtable.RData")


pdf(file="results/ROCNM.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="NM"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))
sigma_values <- sort(sim_table$sigma[scenarios])
for (i in seq_along(sigma_values)){
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(ROC$FPR,
         ROC$TPR,
         type = "l",
         ylab="True positive rate",
         xlab="False positive rate",
         cex.lab=0.75,
         cex.axis=0.6,
         #lty=i,
         col=i)
  }else{
    lines(ROC$FPR,ROC$TPR,
          #lty=i,
          col=i)
  }

}

legend(x=0.6,y=0.5,
       legend=c(expression(sigma *"="* 0.00),
                expression(sigma *"="* 0.50),
                expression(sigma *"="* 0.75),
                expression(sigma *"="* 0.80),
                expression(sigma *"="* 0.85),
                expression(sigma *"="* 0.90),
                expression(sigma *"="* 0.95),
                expression(sigma *"="* 0.99),
                expression(sigma *"="* 1.00)),
       lty=1,
       col=1:9,
       cex=0.75,
       bty="n")

dev.off()




pdf(file="results/ROCSV.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="SV"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))
scenarios <- intersect(scenarios,
                       which(sim_table$initial_frequency==-1))
sigma_values <- sort(sim_table$sigma[scenarios])
for (i in seq_along(sigma_values)){
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(ROC$FPR,
         ROC$TPR,
         type = "l",
         ylab="True positive rate",
         xlab="False positive rate",
         cex.lab=0.75,
         cex.axis=0.6,
         #lty=i,
         col=i)
  }else{
    lines(ROC$FPR,ROC$TPR,
          #lty=i,
          col=i)
  }
  
}

legend(x=0.6,y=0.5,
       legend=c(expression(sigma *"="* 0.00),
                expression(sigma *"="* 0.50),
                expression(sigma *"="* 0.75),
                expression(sigma *"="* 0.80),
                expression(sigma *"="* 0.85),
                expression(sigma *"="* 0.90),
                expression(sigma *"="* 0.95),
                expression(sigma *"="* 0.99),
                expression(sigma *"="* 1.00)),
       lty=1,
       col=1:9,
       cex=0.75,
       bty="n")

dev.off()
