require(dplyr,quietly=T)

# colorblind palette
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")



# load results from running postTreatment/Results.EstimationNe.R script
load(file="results/Ne/Results.RData")

position_x <- c(0.1,0.35,0.65,0.9)
position_y <- c(0.9,0.65,0.35,0.1)
N.values <- sort(unique(results_Ne$N))
(N.values)
sigma.values <- sort(unique(results_Ne$sigma))
(sigma.values)
Ne <- matrix(data=NA,nrow=length(N.values),ncol= length(sigma.values),dimnames=list(N.values,sigma.values))
for (i in seq_along(N.values)){
  Ne[i,] <- (2-sigma.values)*N.values[i]
}
Ne
deltaT.values <- sort(unique(results_Ne$deltaT))
(deltaT.values)

ii <- 3
jj <- 4

old.par <- par(no.readonly=T)
pdf(file="results/FCvsFSTmain.pdf",height=5,width=7)
par(mfrow=c(1,1),mar=c(5.1,5.1,1.1,1.1),oma=c(0,0,0,0))
boxplot( NeFst ~ sigma,
         data = filter(results_Ne,deltaT==deltaT.values[ii],N==N.values[jj]),
         at = sigma.values-0.005,
         outline = FALSE,
         axes = FALSE,
         boxwex = 0.01,
         notch = T,
         col = cbbPalette[5],
         border = cbbPalette[5],
         names = sigma.values,
         ylab = expression("effective population size, "*italic(N[e])),
         xlab = expression("selfing rate, "*sigma),
         ylim = c(0,2000),
         xlim = c(min(sigma.values),max(sigma.values)) )
boxplot( NeFc ~ sigma,
         data = filter(results_Ne,deltaT==deltaT.values[ii],N==N.values[jj]),
         at = sigma.values+0.005,
         outline = FALSE,
         axes = FALSE,
         boxwex = 0.01,
         notch = T,
         col = cbbPalette[7],
         border = cbbPalette[7],
         names = sigma.values,
         ylim = c(0,2000),
         xlim = c(min(sigma.values),max(sigma.values)) ,add=T)
lines(sigma.values,Ne[jj,],col=cbbPalette[1])
axis(side=1,at=sigma.values,las=2,cex.axis=0.8)
axis(2)
box()
legend(x="topleft",
       legend=c(expression(italic(hat(N)[e])*" from "*italic(F)[ST]),
                expression(italic(hat(N)[e])*" from "*italic(F)[C]),
                expression("true "*italic(N[e]))),
       fill=c(cbbPalette[5],
              cbbPalette[7],
              NA),
       border=c(cbbPalette[5],
                cbbPalette[7],
                NA),
       lty=c(NA,NA,1)
       )
dev.off()


pdf(file="results/FCvsFSTsupp.pdf",height=10,width=10)
par(mfrow=c(4,4), mar=c(3,3,1,1), oma=c(6,6,0,0))
first <- T
for (i in seq_along(deltaT.values)){
  for (j in seq_along(N.values)){
    boxplot( NeFst ~ sigma,
             data = filter(results_Ne,deltaT==deltaT.values[i],N==N.values[j]),
             at = sigma.values-0.005,
             outline = FALSE,
             axes = FALSE,
             boxwex = 0.01,
             notch = T,
             col = cbbPalette[5],
             border = cbbPalette[5],
             names = sigma.values,
             ylim = c(0,2000),
             xlim = c(min(sigma.values),max(sigma.values)) )
    boxplot( NeFc ~ sigma,
             data = filter(results_Ne,deltaT==deltaT.values[i],N==N.values[j]),
             at = sigma.values+0.005,
             outline = FALSE,
             axes = FALSE,
             boxwex = 0.01,
             notch = T,
             col = cbbPalette[7],
             border = cbbPalette[7],
             names = sigma.values,
             ylim = c(0,2000),
             xlim = c(min(sigma.values),max(sigma.values)) ,add=T)
    lines(sigma.values,Ne[j,],col=cbbPalette[1])
    axis(side=1,at=sigma.values,las=2,cex.axis=0.5)
    axis(2)
    box()
    if (first==T){
      legend(x="topleft",
             legend=c(expression(italic(hat(N)[e])*" from "*italic(F)[ST]),
                      expression(italic(hat(N)[e])*" from "*italic(F)[C]),
                      expression("true "*italic(N[e]))),
             fill=c(cbbPalette[5],
                    cbbPalette[7],
                    NA),
             border=c(cbbPalette[5],
                      cbbPalette[7],
                      NA),
             lty=c(NA,NA,1))
      first <- F
    }
    
    
    if(i==length(deltaT.values)){
      mtext( substitute(paste(italic(N)*"="*v), list(v=N.values[j])),
             side = 1,
             line = 3,
             cex = 1,
             adj = position_x[j],
             outer = TRUE)
    }
  }
  mtext(substitute(paste(italic(t)*"="*v), list(v=deltaT.values[i])),
        side = 2, 
        line = 3, 
        cex = 1, 
        adj = position_y[i], 
        outer = TRUE)
  
}
mtext(expression("selfing rate, "*sigma), side=1, adj=0.5, cex=1, outer=TRUE)
mtext(expression("effective population size, "*italic(N[e])), side=2, adj=0.5, cex=1, outer=TRUE)
dev.off()

