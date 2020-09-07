# awk '{if ($2=="m1") {print $7, $8}}' populationATequilibrium > ages.txt


ages <- read.table("results/ages.txt")
N <- 500
sigma<-0.95
max_generation <- N*20
ages$V1 <- max_generation-ages$V1
ages$V2 <- ages$V2/N/2

logit <- function(x) log(x/(1-x))

plot(log10(ages$V1),logit(ages$V2))

library(viridis)
library(gplots)
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)

y <- log10(ages$V1+1)
x <- logit(ages$V2)
df <- data.frame(x,y)

expected_age <- function(p){(-2*p)*log(p)/(1-p)}

freqs <- seq(1,2*N-1,1)/N/2
Ne <- (2-sigma)*N/2
times <- expected_age(freqs)*Ne

head(data.frame(log10(times),logit(freqs)))
plot(log10(times),logit(freqs))
pdf(file="results/AlleleAge.pdf",width=9,height=9)
par(mar=c(5.5,5.5,0,0)+0.1)
h2 <- hist2d(df,
             nbins=25,
             col=grey.colors(8,start=0.8,end=0),#col=viridis(10),
             FUN=function(x) log10(length(x)),
             ylim=c(0.0,3.6),
             xlim=c(logit(1/N/2),logit((2*N-1)/N/2)),
             cex.lab=2,
             cex.axis=2,
             ylab=expression(log[10]*"(age)"),
             xlab=expression("logit("*pi*")"))
legend(x=3,y=2,
       legend=c(">1",">3",">10",">30",">100",">300",">1000",">3000"),
       fill=grey.colors(8,start=0.8,end=0),
       cex=2, title="counts",bty="n")

lines(logit(freqs),log10(times),col=7,lwd=4)
dev.off()
h2


