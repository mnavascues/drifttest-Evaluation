source("src/fun/ColorBlindPalette.R")

sigma_values <- c(0.000,0.500,0.750,0.800,0.850,
                  0.900,0.925,0.950,0.975,1.000)



# PLOT FOR SLIM MODEL

pdf(file="results/FCvsFST.pdf",width=4,height=6)
layout(matrix(c(1,2), 2,1,byrow = TRUE))

#par(mfrow=c(2,1),
#    mar=c(2,2,0.2,0.2)+0.1,
#    oma=c(3,3,0,0))

par(mar=c(4,4,0.2,0.2)+0.1)


sim <- 21
simulationID <- "scenario021"
load(file=paste0("results/",simulationID,"/",simulationID,"_results.RData"))
sigma<-0

boxplot( results$NeHatFst/2,
         at = sigma-0.005,
         outline = FALSE,
         boxwex = 0.025,
         notch = T,
         col = cbPalette1[6],
         border = cbPalette1[6],
         names = c("0"),
         ylab = expression("effective population size, "*italic(N)[e]),
         xlab = "",#expression("selfing rate, "*sigma),
         cex.lab=0.75,
         ylim = c(0,1000),
         xlim = c(0,1)  ,cex.axis=0.7,las=1)

boxplot( results$NeHatFc/2,
         at = sigma+0.005,
         axes = FALSE,
         outline = FALSE,
         boxwex = 0.025,
         notch = T,
         col = cbPalette1[7],
         border = cbPalette1[7],
         add=T )





for (sim in 22:30){
  
  if (sim<10)           simulationID <- paste0("scenario00",sim)
  if (sim>=10 & sim<100) simulationID <- paste0("scenario0",sim)
  
  load(file=paste0("results/",simulationID,"/",simulationID,"_results.RData"))
  
  sigma <- sigma_values[sim-20]
  
  boxplot( results$NeHatFst/2,
           at = sigma-0.005,
           outline = FALSE,
           axes = FALSE,
           boxwex = 0.025,
           notch = T,
           col = cbPalette1[6],
           border = cbPalette1[6],
           add=T )
  boxplot( results$NeHatFc/2,
           at = sigma+0.005,
           axes = FALSE,
           outline = FALSE,
           boxwex = 0.025,
           notch = T,
           col = cbPalette1[7],
           border = cbPalette1[7],
           add=T )
  
  
  
}

trueNe <- ((2-sigma_values)*500)/2

sigma_values_axis <- sigma_values
sigma_values_axis[c(7,9)]<-NA

lines(sigma_values,trueNe,lty=2)
axis(side=1,at=sigma_values_axis,las=2,cex.axis=0.6)

text(x=0.95,y=950,label="a",cex=1.5)

box()
#text(0.05,950,"B",cex=1.5)

#mtext(text=expression("selfing rate, "*sigma),
#      side=1,
#      line=1.2,
#      cex=1,
#      outer=TRUE)
#mtext(text=expression("effective population size, "*italic(N)[e]),
#      side=2,
#      line=1.2,
#      cex=1,
#      outer=TRUE)

#legend(x=0.25,y=1000,
#       legend=c(expression(hat(italic(N))[e]*" from "*italic(F)[ST]),
#                expression(hat(italic(N))[e]*" from "*italic(F)[C]),
#                expression("true "*italic(N)[e])),
#       fill=c(cbPalette1[6],
#              cbPalette1[7],
#              NA),
#       border=c(cbPalette1[6],
#                cbPalette1[7],
#                NA),
#       lty=c(NA,NA,2),
#       cex=0.75,
#       bty="n")

#dev.off()




#pdf(file="results/FCvsFSTsuppl.pdf",width=4,height=3)

#par(mfrow=c(2,1),
#    mar=c(2,2,0.2,0.2)+0.1,
#    oma=c(3,3,0,0))

par(mar=c(4,4,0.2,0.2)+0.1)



# PLOT FOR SIMPLE MODEL

sim <- 1
simulationID <- paste0("sim00",sim)
load(file=paste0("results/SimpleModels/",simulationID,".RData"))

boxplot( NeHatFst/2,
         at = sigma-0.005,
         outline = FALSE,
         boxwex = 0.025,
         notch = T,
         col = cbPalette1[6],
         border = cbPalette1[6],
         names = c("0"),
         ylab = expression("effective population size, "*italic(N)[e]),
         xlab = expression("selfing rate, "*sigma),
         cex.lab=0.75,
         ylim = c(0,1000),
         xlim = c(0,1) ,cex.axis=0.7,las=1)

boxplot( NeHatFc/2,
         at = sigma+0.005,
         axes = FALSE,
         outline = FALSE,
         boxwex = 0.025,
         notch = T,
         col = cbPalette1[7],
         border = cbPalette1[7],
         add=T )





for (sim in 2:10){
  
  if (sim<10)           simulationID <- paste0("sim00",sim)
  if (sim>=10 & sim<100) simulationID <- paste0("sim0",sim)
  
  load(file=paste0("results/SimpleModels/",simulationID,".RData"))
  
  boxplot( NeHatFst/2,
           at = sigma,
           outline = FALSE,
           axes = FALSE,
           boxwex = 0.025,
           notch = T,
           col = cbPalette1[6],
           border = cbPalette1[6],
           add=T)
  boxplot( NeHatFc/2,
           at = sigma,
           axes = FALSE,
           outline = FALSE,
           boxwex = 0.025,
           notch = T,
           col = cbPalette1[7],
           border = cbPalette1[7],
           add=T )
  
  
  
}

trueNe <- ((2-sigma_values)*500)/2

sigma_values_axis <- sigma_values
sigma_values_axis[c(7,9)]<-NA

lines(sigma_values,trueNe,lty=2)
axis(side=1,at=sigma_values_axis,las=2,cex.axis=0.6)
box()
#text(0.05,950,"A",cex=1.5)

legend(x=0.25,y=1000,
       legend=c(expression(hat(italic(N))[e]*" from "*italic(F)[ST]),
                expression(hat(italic(N))[e]*" from "*italic(F)[C]),
                expression("true "*italic(N)[e])),
       fill=c(cbPalette1[6],
              cbPalette1[7],
              NA),
       border=c(cbPalette1[6],
                cbPalette1[7],
                NA),
       lty=c(NA,NA,2),
       cex=0.75,
       bty="n")
text(x=0.95,y=950,label="b",cex=1.5)

dev.off()


