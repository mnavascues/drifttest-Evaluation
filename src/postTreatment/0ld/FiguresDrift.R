number_of_loci <- 15
pop_size <- 500
initial_freq  <- rbeta(number_of_loci,0.5,0.5)
initial_count <- rbinom(rep(1,number_of_loci),pop_size,initial_freq)

dT <-25

trajectories <- matrix(NA,dT+1,number_of_loci)
trajectories[1,] <- initial_count

for (generation in 1:dT+1){
  
  freq <- trajectories[generation-1,]/pop_size
  
  trajectories[generation,] <- rbinom(rep(1,number_of_loci),pop_size,freq)
  
}
trajectories[1,1]<-50
for (generation in 1:dT+1){
  freq <- trajectories[generation-1,1]/pop_size
  if (trajectories[generation-1,1]/pop_size==1){
    trajectories[generation,1] <- pop_size
  }else{
    trajectories[generation,1] <- rbinom(1,pop_size,freq*1.1)
    if (trajectories[generation,1]>=pop_size) {
      trajectories[generation,1] <- pop_size
    }
  }
}





plot(c(0,dT),
     c(trajectories[1,1]/pop_size,trajectories[dT+1,1]/pop_size),
     type="p",
     ylim=c(0,1),
     ylab="Allele frequency",xlab="time",
     pch=1)
for (locus in 2:number_of_loci){
  points(0,trajectories[1,locus]/pop_size,pch=locus)
  points(dT,trajectories[dT+1,locus]/pop_size,pch=locus)
}
points(0,trajectories[1,1]/pop_size,col="red",pch=1,lwd=3)
points(dT,trajectories[dT+1,1]/pop_size,col="red",pch=1,lwd=3)





plot(0:25,trajectories[,1]/pop_size,
     type="l",
     ylim=c(0,1),
     ylab="Allele frequency",xlab="time",
     col="grey",lwd=3)
points(0,trajectories[1,1]/pop_size,pch=1)
points(dT,trajectories[dT+1,1]/pop_size,pch=1)
for (locus in 2:number_of_loci){
  lines(0:25,trajectories[,locus]/pop_size,col="grey",lwd=2)
  points(0,trajectories[1,locus]/pop_size,pch=locus)
  points(dT,trajectories[dT+1,locus]/pop_size,pch=locus)
}
lines(0:25,trajectories[,1]/pop_size,col="red",lwd=2)
points(0,trajectories[1,1]/pop_size,col="red",pch=1,lwd=3)
points(dT,trajectories[dT+1,1]/pop_size,col="red",pch=1,lwd=3)

number_of_simulations <- 1000

simulated_trajectories <- matrix(NA,dT+1,number_of_simulations)

simulated_trajectories[1,] <- trajectories[1,1]

for (generation in 1:dT+1){
  
  freq <- simulated_trajectories[generation-1,]/pop_size
  
  simulated_trajectories[generation,] <- rbinom(rep(1,number_of_simulations),pop_size,freq)
  
}


plot(0:dT,simulated_trajectories[,1]/pop_size,
     type="l",
     ylim=c(0,1),
     ylab="Allele frequency",xlab="time",
     col="grey",lwd=1)
points(dT,simulated_trajectories[dT+1,1]/pop_size)
for (locus in 2:number_of_simulations){
  lines(0:dT,simulated_trajectories[,locus]/pop_size,col="grey",lwd=1)
  points(dT,simulated_trajectories[dT+1,locus]/pop_size)
  
}
points(0,simulated_trajectories[1,1]/pop_size)

#lines(0:dT,trajectories[,1]/pop_size, col="orange",lwd=3)
points(0,trajectories[1,1]/pop_size,col="red")
points(dT,trajectories[dT+1,1]/pop_size,col="red")

