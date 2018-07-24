source("src/fun/F_stats_tools.R")
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)

Ne     <- 500
dT <- 25
sample_size <- 100

initial.allele.freq <- 0.05
final.allele.freq   <- 0.20

observed_counts <- rbind(c(sample_size*initial.allele.freq,
                           sample_size*(1-initial.allele.freq),
                           sample_size*final.allele.freq,
                           sample_size*(1-final.allele.freq)))
observed_Fst <- single_locus_compute_FST_from_counts(observed_counts)

num_of_rep<-2000
sim_Fst <- array(NA,num_of_rep)
for (rep in 1:num_of_rep){
  trajectory <- array(NA,dT+1)
  trajectory[1] <- rbeta( n=1, shape1=1+sample_size*initial.allele.freq,
                               shape2=1+sample_size*(1-initial.allele.freq) )
  for (g in 1:dT) {
    trajectory[g+1] <- (rbinom(1,Ne,trajectory[g])) / Ne
  }

  if (rep==1){
    pdf(file="results/Test.pdf",width=6, height=4.5)
    plot(0:dT,
         trajectory,
         xlim=c(0,dT),
         ylim=c(0,0.5),
         ylab="Allele frequency",xlab="Time",col=rgb(0,0,0,0.1),type="l")
  }else{
    if (trajectory[dT]>0.01){
      if (rep<1000) lines(0:dT,trajectory,col=rgb(0,0,0,0.1))
      count_t1 <- rbinom(1,sample_size,trajectory[1])
      count_t2 <- rbinom(1,sample_size,trajectory[dT])
      counts <- rbind(c(count_t1,sample_size-count_t1,
                        count_t2,sample_size-count_t2))
      sim_Fst[rep] <- single_locus_compute_FST_from_counts(counts)
    }
  }
}
points(x=c(0,dT),
       y=c(initial.allele.freq,y1=final.allele.freq),
       col=2,cex=1.8,pch=16)
dev.off()

pdf(file="results/NullDistribution.pdf",width=6, height=4.5)
plot(density(sim_Fst,na.rm=T),
     xlim=c(-0.03,0.2),
     ylab="Probability density",
     xlab=expression(italic(F)[ST]),
     main="")
abline(v=observed_Fst,col=2)

x <- density(sim_Fst,na.rm=T)$x[which(density(sim_Fst,na.rm=T)$x>observed_Fst)]
y <- density(sim_Fst,na.rm=T)$y[which(density(sim_Fst,na.rm=T)$x>observed_Fst)]

p_value <- sum(sim_Fst>observed_Fst,na.rm=T)/sum(! is.na(sim_Fst))

text(0.1,25,expression("observed "*italic(F)[ST]*"=0.09"),pos=4,cex=1.2,col=2)
text(0.1,20,expression(italic(p)*"-value=0.03"),pos=4,cex=1.2,col=3)


polygon(c(observed_Fst,x,30),c(0,y,0),col=3,border=3)
dev.off()

