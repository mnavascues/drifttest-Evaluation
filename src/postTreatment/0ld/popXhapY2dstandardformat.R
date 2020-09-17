source("F_stats_tools.R")


# read simulations description table
sim_table            <- read.table("Simulations1locus/simparams.txt",header=T)
number_of_replicates <- 100000
N                    <- sim_table$population_size                 # number of individuals
dT                   <- sim_table$selection_period_duration
scenarios            <- as.character(sim_table$simID)
selfing_coeff        <- sim_table$sigma 
Ne                   <- (2-selfing_coeff)*N  # number of gene copies
Exp_Fst              <- dT/(dT+2*Ne)
E_Fc                 <- 1/100 + 1/100 + dT/Ne - 2/Ne

(sim_table)
unique(N);unique(dT);unique(scenarios);unique(selfing_coeff);unique(Ne);unique(Exp_Fst);unique(E_Fc)

num_of_loci <- 1000
num_of_lots <- round(number_of_replicates/num_of_loci)

estimated_Fst_1 <- matrix(NA,nrow=num_of_lots,ncol=length(scenarios))

for (scen in seq_along(scenarios)){
  replic<-1
  for (lote in seq_len(num_of_lots)){
    genotype <- matrix(c(1:100,rep(1,50),rep(2,50)),nrow=100,ncol=2)
    
    for (locus in seq_len(num_of_loci)){
      success <- try( load(paste0("Simulations1locus/",scenarios[scen],"/",scenarios[scen],"_",replic,"_data.RData")) , silent=T)
      if(class(success)=="try-error"){
        cat(paste0("Simulation ",scenarios[scen]," replicate ",replic," file does not exist\n"))
      }else{
        cat(paste0("Simulation ",scenarios[scen]," replicate ",replic,"\n"))
        genotype <- cbind(genotype,SNP_data$genotype_data[,3])
      }
      replic<-replic+1  
    }
    if(ncol(genotype)>1000){
      estimated_Fst_1[lote,scen] <- compute_Fstats(genotype)$F_ST 
    }
    
    #compute_Fstats(genotype)$freq
  }
}


boxplot(estimated_Fst_1[sample(1:100,20),],
        names=selfing_coeff,
        ylim=c(0.01,0.035),
        ylab=expression(italic("F")[ST]),
        xlab="Selfing Rate",
        main="SLiM")
lines(seq_along(scenarios), Exp_Fst, col="red")
points(seq_along(scenarios), Exp_Fst, col="red")
points(seq_along(scenarios), apply(estimated_Fst_1,2,mean,na.rm=T), col="blue")

save(estimated_Fst_1,file="Simulations1locus/estimated_Fst.RData")

#boxplot(estimated_Fst_1[,c(1:3,6,8:10)],names=selfing_coeff[c(1:3,6,8:10)],ylim=c(0.01,0.030))
#lines(seq_along(scenarios[c(1:3,6,8:10)]), Exp_Fst[c(1:3,6,8:10)], col="red")
#points(seq_along(scenarios[c(1:3,6,8:10)]), Exp_Fst[c(1:3,6,8:10)], col="red")
#points(seq_along(scenarios[c(1:3,6,8:10)]), apply(estimated_Fst_1[,c(1:3,6,8:10)],2,mean), col="blue")















# read simulations description table
sim_table            <- read.table("SimulationsFst/simparams.txt",header=T)
number_of_replicates <- 200
N                    <- sim_table$population_size                 # number of individuals
dT                   <- sim_table$selection_period_duration
scenarios            <- as.character(sim_table$simID)
selfing_coeff        <- sim_table$sigma 
Ne                   <- (2-selfing_coeff)*N  # number of gene copies
Exp_Fst              <- dT/(dT+2*Ne)
E_Fc                 <- 1/100 + 1/100 + dT/Ne - 2/Ne

(sim_table)
unique(N);unique(dT);unique(scenarios);unique(selfing_coeff);unique(Ne);unique(Exp_Fst);unique(E_Fc)


estimated_Fst_2 <- matrix(NA,nrow=number_of_replicates,ncol=length(scenarios))

for (scen in seq_along(scenarios)){
  for (replic in seq_len(number_of_replicates)){

      success <- try( load(paste0("SimulationsFst/",scenarios[scen],"/",scenarios[scen],"_",replic,"_F_stats.RData")) , silent=T)
      if(class(success)=="try-error"){
        cat(paste0("Simulation replicate ",replic," file does not exist\n"))
      }else{
        cat(paste0("Simulation replicate ",replic,"\n"))
        estimated_Fst_2[replic,scen] <- Fst_hat
      }
  }
}
boxplot(estimated_Fst_2,names=selfing_coeff,ylim=c(-0.02,0.35))
lines(seq_along(scenarios), Exp_Fst, col="red")
points(seq_along(scenarios), Exp_Fst, col="red")
points(seq_along(scenarios), apply(estimated_Fst_2,2,mean), col="blue")


