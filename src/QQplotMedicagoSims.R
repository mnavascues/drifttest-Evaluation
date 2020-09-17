arabidopsis_p_value <- read.table("data/arabidopsis.txt", header=T)
head(arabidopsis_p_value)
medicago_p_value <- readRDS(file="results/medicago_p_value.rds")
load(file="results/scenario002/scenario002_results.RData")
res_sims_100 <- results_per_locus
load(file="results/scenario009/scenario009_results.RData")
res_sims_99 <- results_per_locus
load(file="results/scenario007/scenario007_results.RData")
res_sims_95 <- results_per_locus
load(file="results/scenario003/scenario003_results.RData")
res_sims_90 <- results_per_locus
load(file="results/scenario004/scenario004_results.RData")
res_sims_80 <- results_per_locus

source("src/fun/ColorBlindPalette.R")



library("qqman")
citation("qqman")



qq2 <- function (pvector, ...) 
{
  if (!is.numeric(pvector)) 
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                       is.finite(pvector) & pvector < 1 & pvector > 0]
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  def_args <- list(pch = 20, xlim = c(0, max(e)), ylim = c(0, 
                                                           max(o)), xlab = expression(Expected ~ ~-log[10](italic(p))), 
                   ylab = expression(Observed ~ ~-log[10](italic(p))))
  dotargs <- list(...)
  tryCatch(do.call("plot", c(list(x = e, y = o), def_args[!names(def_args) %in% 
                                                            names(dotargs)], dotargs)), warn = stop)
  abline(0, 1, col = "black", lty=2)
}


pdf(file="results/qqplot_Medicago.pdf",width=4.5,height=9)
par(mar=c(4,4,0.2,1)+0.1)
#layout(matrix(1:6, 2,3,byrow = TRUE))
layout(matrix(1:2, 2,1,byrow = TRUE))
qq2(medicago_p_value,xlim=c(0,6),ylim=c(0,6),type="l",col=cbPalette1[2],lwd=2)

pvector<-arabidopsis_p_value$p_value
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                     is.finite(pvector) & pvector < 1 & pvector > 0]
o <- -log10(sort(pvector, decreasing = FALSE))
e <- -log10(ppoints(length(pvector)))
lines(e,o,col=cbPalette1[3],lwd=2)

legend(0,6,
       legend=c(expression(italic(Medicago)*" "*italic(truncatula)),
                expression(italic(Arabidopsis)*" "*italic(thaliana))),
       lwd=2,bty="n",
       col=cbPalette1[2:3])

text(5.8,5.8,"a",cex=2)

qq2(res_sims_100$p_value,xlim=c(0,6),ylim=c(0,6),type="l",col=cbPalette1[4],lwd=2)

pvector<-res_sims_99$p_value
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                     is.finite(pvector) & pvector < 1 & pvector > 0]
o <- -log10(sort(pvector, decreasing = FALSE))
e <- -log10(ppoints(length(pvector)))
lines(e,o,col=cbPalette1[5],lwd=2)

pvector<-res_sims_95$p_value
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                     is.finite(pvector) & pvector < 1 & pvector > 0]
o <- -log10(sort(pvector, decreasing = FALSE))
e <- -log10(ppoints(length(pvector)))
lines(e,o,col=cbPalette1[6],lwd=2)

pvector<-res_sims_90$p_value
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                     is.finite(pvector) & pvector < 1 & pvector > 0]
o <- -log10(sort(pvector, decreasing = FALSE))
e <- -log10(ppoints(length(pvector)))
lines(e,o,col=cbPalette1[7],lwd=2)

pvector<-res_sims_80$p_value
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                     is.finite(pvector) & pvector < 1 & pvector > 0]
o <- -log10(sort(pvector, decreasing = FALSE))
e <- -log10(ppoints(length(pvector)))
lines(e,o,col=cbPalette1[8],lwd=2)


legend(0,6,
       legend=c(expression(sigma*"=1.00"),
                expression(sigma*"=0.99"),
                expression(sigma*"=0.95"),
                expression(sigma*"=0.90"),
                expression(sigma*"=0.80")),
       lwd=2,bty="n",
       col=cbPalette1[4:8],
       title="Simulations")


text(5.8,5.8,"b",cex=2)

dev.off()


  

