
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)

load(file="results/simtable.RData")


pdf(file="results/ROC.pdf",width=4,height=8)
layout(matrix(c(1,2), 2,1,byrow = TRUE))

par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="NM"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))
sigma_values <- sim_table$sigma[scenarios]

scenarios <- c(1,2,3,4,6,7)
sigma_values <- sort(sim_table$sigma[scenarios])

for (i in seq_along(sigma_values)){
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  neutral_loci <- intersect(which(results_per_locus$neutral),
                            which(!is.na(results_per_locus$p_value)))
  selected_loci <- intersect(intersect(which(!results_per_locus$neutral),
                                       which(results_per_locus$distance==0)), #results_per_locus$centimorgan<=1
                             which(!is.na(results_per_locus$p_value)))
  threshold4ROC <- c(0,10^-seq(10,3,-0.1),seq(0.002,1,0.001))
  FPR <- TPR <- array(NA,length(threshold4ROC))
  for (T.hold in seq_along(threshold4ROC)){
    FPR[T.hold] <- sum(results_per_locus$p_value[neutral_loci] <= threshold4ROC[T.hold]) / length(neutral_loci)
    TPR[T.hold] <- sum(results_per_locus$p_value[selected_loci] <= threshold4ROC[T.hold]) / length(selected_loci)
  }
  ROC<-data.frame(FPR,TPR)
  #ROC <- rbind(c(0,0),ROC,c(1,1))
  
  #head(ROC)
  
  if (i==1){
    plot(ROC$FPR,
         ROC$TPR,
         type = "l",
         ylab="True positive rate",
         xlab="",
         xlim=c(0,1),
         ylim=c(0,1),
         cex.lab=0.75,
         cex.axis=0.6,
         #lty=i,
         col=i)
  }else{
    lines(ROC$FPR,ROC$TPR,
          #lty=i,
          col=i)
  }

  points(FPR_thresholdPvalue0.001,POWER_thresholdPvalue0.001,col=i)
  
  
  
}

legend(x=0.6,y=0.5,
       legend=c(expression(sigma *"="* 0.00),
                expression(sigma *"="* 0.50),
                expression(sigma *"="* 0.80),
                expression(sigma *"="* 0.90),
                expression(sigma *"="* 0.95),
                expression(sigma *"="* 1.00)),
       lty=1,
       col=1:6,
       cex=0.75,
       bty="n")
text(x=1,y=1,label="a",cex=1.5)

#dev.off()




#pdf(file="results/ROCSV.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="SV"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))
scenarios <- intersect(scenarios,
                       which(sim_table$initial_frequency==-1))
sigma_values <- sim_table$sigma[scenarios]
scenarios<-c(10,11,13,15,16,18)
sigma_values <- sort(sim_table$sigma[scenarios])
for (i in seq_along(sigma_values)){
  
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  neutral_loci <- intersect(which(results_per_locus$neutral),
                            which(!is.na(results_per_locus$p_value)))
  selected_loci <- intersect(intersect(which(!results_per_locus$neutral),
                                       which(results_per_locus$distance==0)), #results_per_locus$centimorgan<=1
                             which(!is.na(results_per_locus$p_value)))
  threshold4ROC <- c(0,10^-seq(10,3,-0.1),seq(0.002,1,0.001))
  FPR <- TPR <- array(NA,length(threshold4ROC))
  for (T.hold in seq_along(threshold4ROC)){
    FPR[T.hold] <- sum(results_per_locus$p_value[neutral_loci] <= threshold4ROC[T.hold]) / length(neutral_loci)
    TPR[T.hold] <- sum(results_per_locus$p_value[selected_loci] <= threshold4ROC[T.hold]) / length(selected_loci)
  }
  ROC<-data.frame(FPR,TPR)
  
  if (i==1){
    plot(ROC$FPR,
         ROC$TPR,
         type = "l",
         ylab="True positive rate",
         xlab="False positive rate",
         xlim=c(0,1),
         ylim=c(0,1),
         cex.lab=0.75,
         cex.axis=0.6,
         #lty=i,
         col=i)
  }else{
    lines(ROC$FPR,ROC$TPR,
          #lty=i,
          col=i)
  }
  points(FPR_thresholdPvalue0.001,POWER_thresholdPvalue0.001,col=i)
  
}

#legend(x=0.6,y=0.5,
#       legend=c(expression(sigma *"="* 0.00),
#                expression(sigma *"="* 0.50),
#                expression(sigma *"="* 0.80),
#                expression(sigma *"="* 0.90),
#                expression(sigma *"="* 0.95),
#                expression(sigma *"="* 1.00)),
#       lty=1,
#       col=1:6,
#       cex=0.75,
#       bty="n")
text(x=1,y=1,label="b",cex=1.5)

dev.off()
