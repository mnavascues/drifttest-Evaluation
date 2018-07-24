
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)

#r  <- 1e-8 
#kf <- 1/2/r
#d <- 50*log(1/(1-2*4e7*r))




load(file="results/simtable.RData")


pdf(file="results/SelectionFootprintNM.pdf",width=4,height=3)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="NM"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))

scenarios <- c(1,6,3,9)
sigma_values <- sort(sim_table$sigma[scenarios])

for (i in seq_along(sigma_values)){
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(SelectionFootprint_thresholdPvalue0.001[,"distance"],
         SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,1),
         #xlim=c(0,40),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)")
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)

  }else{
    lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
          SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
          #lty=i,
          col=i)
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }
  
}

legend(x=1,y=1,
       legend=c(expression(sigma *"="* 0.00),
                expression(sigma *"="* 0.50),
#                expression(sigma *"="* 0.75),
#                expression(sigma *"="* 0.80),
#                expression(sigma *"="* 0.85),
                expression(sigma *"="* 0.90),
#                expression(sigma *"="* 0.95),
#                expression(sigma *"="* 0.99),
                expression(sigma *"="* 0.99)),
       lty=1,
       col=1:4,
       cex=0.75,
       bty="n")

dev.off()





pdf(file="results/SelectionFootprintSelCoef.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(49,50,51,52,1)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(SelectionFootprint_thresholdPvalue0.001[,"distance"],
         SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,1),
         #xlim=c(0,40),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)")
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }else{
    lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
          SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
          #lty=i,
          col=i)
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }
  
}

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






























pdf(file="results/SelectionFootprintTime1.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(49,53,54,59)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(SelectionFootprint_thresholdPvalue0.001[,"distance"],
         SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,1),
         #xlim=c(0,40),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)")
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }else{
    lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
          SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
          #lty=i,
          col=i)
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }
  
}

legend(x=1,y=1,
       legend=c(expression(tau *"="* 25),
                expression(tau *"="* 50),
                expression(tau *"="* 100),
                expression(tau *"="* 200)),
       lty=1,
       col=1:4,
       cex=0.75,
       bty="n")

dev.off()













pdf(file="results/SelectionFootprintTime5.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(58,57,1,55,56,60)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(SelectionFootprint_thresholdPvalue0.001[,"distance"],
         SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,1),
         #xlim=c(0,40),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)")
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }else{
    lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
          SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
          #lty=i,
          col=i)
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }
  
}

legend(x=1,y=1,
       legend=c(expression(tau *"="* 5),
                expression(tau *"="* 10),
                expression(tau *"="* 25),
                expression(tau *"="* 50),
                expression(tau *"="* 100),
                expression(tau *"="* 200)),
       lty=1,
       col=1:6,
       cex=0.75,
       bty="n")

dev.off()











palette(cbPalette1)

pdf(file="results/SelectionFootprintSV.pdf",width=4,height=3)
par(mar=c(4,4,0.2,0.2)+0.1)


#plot new mutation (scenario007)
load(file=paste0("results/scenario007/scenario007_results.RData"))
plot(SelectionFootprint_thresholdPvalue0.001[,"distance"],
     SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
     cex.lab=0.75,
     cex.axis=0.6,
     col="white",#1,
     ylim=c(0,1),
     #xlim=c(0,40),
     log="x",
     type="l",
     ylab="Proportion of positives",
     xlab="distance to selected locus (cM)")
#points(60,FPR_thresholdPvalue0.001,col=1)
#points(0.009,POWER_thresholdPvalue0.001,col=1)



scenarios <- c(31,35,39)
scenarios <- c(31,39)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
          SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
          #lty=i,
          col=i+5)
  points(60,FPR_thresholdPvalue0.001,col=i+5)
  points(0.009,POWER_thresholdPvalue0.001,col=i+5)
  
    
}
scenarios <- c(40,44,48)
scenarios <- c(40,48)
for (i in seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
        SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
        lty=2,
        col=i+5)
  points(60,FPR_thresholdPvalue0.001,col=i+5,pch=4)
  points(0.009,POWER_thresholdPvalue0.001,col=i+5,pch=4)
}





legend(x=0.05,y=1,
       legend=c(#expression("new mutation"),
                expression("ancestral allele, "*pi[0]*"="* 0.1),
                expression("derived allele, "*pi[0]*"="* 0.1),
                #expression("ancestral allele, "*pi[0]*"="* 0.5),
                #expression("derived allele, "*pi[0]*"="* 0.5),
                expression("ancestral allele, "*pi[0]*"="* 0.9),
                expression("derived allele, "*pi[0]*"="* 0.9)),
       lty=c(1,2,1,2),#c(1,1,2,1,2,1,2),
       pch=c(1,4,1,4),#c(1,1,4,1,4,1,4), 
       col=c(6,6,7,7),#c(1,2,2,3,3,4,4),
       cex=0.75,
       bty="n")

dev.off()

