source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)


library(qqman)
load(file="results/simtable.RData")




load(file=paste0("results/scenario021/scenario021_results.RData"))
p_valueMAF <- results_per_locus$p_value[!is.na(results_per_locus$p_value)]
uniform.quantilesMAF <- qunif((1:length(p_valueMAF))/(length(p_valueMAF) +1))

load(file=paste0("results/scenario021b/scenario021b_results.RData"))
p_value <- results_per_locus$p_value[!is.na(results_per_locus$p_value)]
uniform.quantiles <- qunif((1:length(p_value))/(length(p_value) +1))





pdf(file="results/QQplotMAF.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)


# CI: http://www.gettinggeneticsdone.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html

plot(c(-1,2),
     c(-1,2),
     cex.lab=0.75,
     cex.axis=0.6,
     type = "l",
     lty=2,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = expression(Expected ~ ~(italic(p))),
     ylab = expression(Observed ~ ~(italic(p))),
     col=1,
     main = "")

lines(sort(uniform.quantilesMAF),
     sort(p_valueMAF),col=6,lwd=2)
lines(sort(uniform.quantiles),
      sort(p_value),col=7,lwd=2)


par(new = TRUE)
split.screen(matrix(data= c(grconvertX(c(0,0.5), from="user", to="ndc"),
                            grconvertY(c(0.6,1), from="user", to="ndc"),
                            grconvertX(c(0.5,1), from="user", to="ndc"),
                            grconvertY(c(0,0.4), from="user", to="ndc")),
                    ncol=4,byrow = T))
screen(1)
par(mar=c(1,0,0,0))
hist(p_value,main="",xlab="",ylab="",freq=F, col=7,border=F,axes=F)
axis(side=1,tck=-0.04,labels=F)
axis(side=1,cex.axis=0.5,tick=F,lwd=-1,line=-1)
mtext("p-value", side=1, line=0.5, cex=0.5, adj=0.5)
box()
screen(2)
par(mar=c(1,0,0,0))
hist(p_valueMAF,main="",xlab="",ylab="",freq=F, col=6,border=F,axes=F)
axis(side=1,tck=-0.04,labels=F)
axis(side=1,cex.axis=0.5,tick=F,lwd=-1,line=-1)
mtext("p-value", side=1, line=0.5, cex=0.5, adj=0.5)
box()
dev.off()

