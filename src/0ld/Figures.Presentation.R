num_of_loci <- 15
initial.allele.freq <- rbeta(num_of_loci,0.9,0.9)
Ne     <- 200
dT <- 20
final.allele.freq <- array(NA, num_of_loci)
initial.allele.freq[1] <- 0.99
initial.allele.freq[2] <- 0.00
final.allele.freq[1]   <- 1.00
final.allele.freq[2]   <- 0.01
for (locus in 3:num_of_loci){
  allele.freq <- initial.allele.freq[locus]
  for (g in 1:dT) {
    allele.freq <- (rbinom(1,Ne,allele.freq)) / Ne
  }
  final.allele.freq[locus] <- allele.freq  
  
}
(final.allele.freq)

initial.allele.freq <- c(initial.allele.freq,0.05)
final.allele.freq   <- c(final.allele.freq,0.60)

pdf(file="Data.pdf",width=6, height=4.5)
plot(c(0,dT),
     c(initial.allele.freq[num_of_loci+1],final.allele.freq[num_of_loci+1]),
     type="n",
     xlim=c(0,20),
     ylim=c(0,1),
     ylab="Allele frequency",xlab="Time")
for(locus in seq_len(num_of_loci+1)){
  arrows(x0=0,x1=dT,
         y0=initial.allele.freq[locus],y1=final.allele.freq[locus],
         length=0.1)
}

dev.off()

pdf(file="DataHighlight.pdf",width=6, height=4.5)
plot(c(0,dT),
     c(initial.allele.freq[num_of_loci+1],final.allele.freq[num_of_loci+1]),
     type="n",
     xlim=c(0,20),
     ylim=c(0,1),
     ylab="Allele frequency",xlab="Time",pch=12,col="red")
for(locus in seq_len(num_of_loci)){
  arrows(x0=0,x1=dT,
         y0=initial.allele.freq[locus],y1=final.allele.freq[locus],
         length=0.1,col="grey")
}
arrows(x0=0,x1=dT,
       y0=initial.allele.freq[num_of_loci+1],y1=final.allele.freq[num_of_loci+1],
       length=0.1,col="red")
dev.off()


pdf(file="DataHighlight2.pdf",width=6, height=4.5)
plot(c(0,dT),
     c(initial.allele.freq[num_of_loci+1],final.allele.freq[num_of_loci+1]),
     type="n",
     xlim=c(0,20),
     ylim=c(0,1),
     ylab="Allele frequency",xlab="Time")
for(locus in 1:2){
  arrows(x0=0,x1=dT,
         y0=initial.allele.freq[locus],y1=final.allele.freq[locus],
         length=0.1,col="blue")
}
for(locus in 3:(num_of_loci+1)){
  arrows(x0=0,x1=dT,
         y0=initial.allele.freq[locus],y1=final.allele.freq[locus],
         length=0.1)
}
dev.off()


#test OLD

num_of_rep<-1000
for (rep in 1:num_of_rep){
  trajectory <- array(NA,dT+1)
  trajectory[1] <- 0.05
  for (g in 1:dT) {
    trajectory[g+1] <- (rbinom(1,Ne,trajectory[g])) / Ne
  }
  
  if (rep==1){
    pdf(file="TestOld.pdf",width=6, height=4.5)
    plot(0:dT,
         trajectory,
         xlim=c(0,20),
         ylim=c(0,1),
         ylab="Allele frequency",xlab="Time",col=rgb(0,0,0,0.1),type="l")
    
  }else{
    lines(0:dT,trajectory,col=rgb(0,0,0,0.1))
  }
  
  
}
arrows(x0=0,x1=dT,
       y0=initial.allele.freq[num_of_loci+1],y1=final.allele.freq[num_of_loci+1],
       length=0.1,col="red")
dev.off()






#test NEW

num_of_rep<-1000
for (rep in 1:num_of_rep){
  trajectory <- array(NA,dT+1)
  trajectory[1] <- rbeta( n=1, shape1=1+5, shape2=1+95 )
  for (g in 1:dT) {
    trajectory[g+1] <- (rbinom(1,Ne,trajectory[g])) / Ne
  }

  if (rep==1){
    pdf(file="TestNew.pdf",width=6, height=4.5)
    plot(0:dT,
         trajectory,
         xlim=c(0,20),
         ylim=c(0,1),
         ylab="Allele frequency",xlab="Time",col=rgb(0,0,0,0.1),type="l")
    pdf(file="TestNewMAF.pdf",width=6, height=4.5)
    plot(0:dT,
         trajectory,
         xlim=c(0,20),
         ylim=c(0,1),
         ylab="Allele frequency",xlab="Time",col=rgb(0,0,0,0.1),type="l")
    dev.set(dev.prev())
  }else{
    lines(0:dT,trajectory,col=rgb(0,0,0,0.1))
    dev.set(dev.next())
    if (trajectory[dT]>0.01) lines(0:dT,trajectory,col=rgb(0,0,0,0.1))
    dev.set(dev.prev())
  }
  
    
}
arrows(x0=0,x1=dT,
       y0=initial.allele.freq[num_of_loci+1],y1=final.allele.freq[num_of_loci+1],
       length=0.1,col="red")
dev.set(dev.next())
arrows(x0=0,x1=dT,
       y0=initial.allele.freq[num_of_loci+1],y1=final.allele.freq[num_of_loci+1],
       length=0.1,col="red")
dev.off()
dev.off()



