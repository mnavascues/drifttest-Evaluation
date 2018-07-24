
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)


bp2centimorgan <- function(distance_in_bp,recombination_rate){
  d <- 50*log(1/(1-2*distance_in_bp*recombination_rate))
  return(d)
}
r<-1e-8

sigma <- c(0,0.5,0.9,0.99)
dT    <- 25
N     <- 500
Ne    <- (2-sigma)*N
e_Fst <- dT/(dT+4*Ne) 





load(file="results/simtable.RData")

pdf(file="results/FstDistance.pdf",width=6,height=4)
par(mar=c(4,4,0.2,1)+0.1)

load(file=paste0("results/scenario001/scenario001_results.RData"))
boxplot(Fst ~ distance,
        data=fst_with_distance,
        axes=F,
        ylim=c(0,0.6),
        xlim=c(0,40),
        notch = T,
        boxwex = 0.4,
        xaxp = c(0,10,9),
        outline = FALSE,
        col = rgb(1,1,1,max=255,alpha=15*255/100),
        border = rgb(1,1,1,max=255,alpha=15*255/100),
        #cex.axis=0.7,
        ylab = expression(italic(F)[ST]),
        staplewex=0,
        las=1,
        xlab = expression("distance to selected site (bp)"))

tick_names <- unique(sort(fst_with_distance$distance))
tick_names <- tick_names[which(tick_names%%10000000==0)]
#tick_names <- round(bp2centimorgan(tick_names,r))
axis(side=1,
     at=seq(1,101,10),
     labels = tick_names,
     cex.axis=1)
axis(side=2,
     at=seq(0,0.8,0.2),
     cex.axis=1)
box()

means <- tapply(fst_with_distance$Fst,fst_with_distance$distance,mean,na.rm=T)
lines(means,col="black",lwd=2)
abline(h=e_Fst[1],col="black",lty=2)


scenarios <- c(6,3,9)
for (i in seq_along(scenarios)){
     
  simID <- paste0("scenario00",scenarios[i])
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))

  transparent <- rgb( col2rgb(cbPalette1[i+1])[1],
                      col2rgb(cbPalette1[i+1])[2],
                      col2rgb(cbPalette1[i+1])[3],
                      max=255,
                      alpha = 15*255/100)
  
  boxplot(fst_with_distance$Fst ~ fst_with_distance$distance,
          ylim=c(0,1),
          notch = T,
          boxwex = 0.4,
          outline = FALSE,
          staplewex=0,
          axes = FALSE,
          col = transparent,
          border = transparent,
          add=T)
  means <- tapply(fst_with_distance$Fst,fst_with_distance$distance,mean,na.rm=T)
  lines(means,col=cbPalette1[i+1],lwd=2)
  abline(h=e_Fst[i+1],col=cbPalette1[i+1],lty=2)
  
}

dev.off()

