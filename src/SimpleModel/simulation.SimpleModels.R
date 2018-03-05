# Simulation under simple models
#################################

number_of_replicates <- 100 
number_of_loci       <- 10000
N                    <- 500
sample_size          <- c(50,50)
dT                   <- 25
seed4random          <- 166656
set.seed(seed4random)

source("src/fun/F_stats_tools.R")

sigma_values <- c(0,0.5,0.75,0.8,0.85,0.9,0.925,0.95,0.975,1)

for (sim in 1:10){
  simulationID <- NA
  if (sim<10)            simulationID <- paste0("sim00",sim)
  if (sim>=10 & sim<100) simulationID <- paste0("sim0",sim)
  
  
  model     <- "Constant" # "Constant" # "Beta"
  parameter <- 0.5
  sigma     <- sigma_values[sim]
  Fis       <- sigma/(2-sigma)
  Ne        <- round((2-sigma)*N)
  
  NeHatFc        <- NeHatFst        <- array(NA,number_of_replicates)
  NeHatFcMAF0.1  <- NeHatFstMAF0.1  <- array(NA,number_of_replicates)
  NeHatFcMAF0.05 <- NeHatFstMAF0.05 <- array(NA,number_of_replicates)
  NeHatFcMAF0.01 <- NeHatFstMAF0.01 <- array(NA,number_of_replicates)
  
  if (Fis==0){
    for ( r in seq_len(number_of_replicates) ){
      if (r==1)     print(paste(simulationID,"Replicate",r))
      if (r%%50==0) print(paste(simulationID,"Replicate",r))
      
      if (model=="Constant") px <- rep(parameter,number_of_loci)   
      if (model=="Beta")     px <- rbeta(number_of_loci,parameter,parameter)    
      py <- rep(NA,number_of_loci)
      
      genotype_data <- cbind( 1:sum(sample_size),
                              c(rep(1,sample_size[1]),rep(2,sample_size[1])) )
      
      for (locus in seq_len(number_of_loci) ){
        f.drift <- px[locus]
        for (g in 1:dT) {
          f.drift <- (rbinom(1,Ne,f.drift)) / Ne
        }
        py[locus] <- f.drift
        
        geno <- c(colSums( rbind(sample(c(0,1),sample_size[1],replace=TRUE,prob=c(px[locus],1-px[locus])),
                                 sample(c(0,1),sample_size[1],replace=TRUE,prob=c(px[locus],1-px[locus])))),
                  colSums( rbind(sample(c(0,1),sample_size[1],replace=TRUE,prob=c(py[locus],1-py[locus])),
                                 sample(c(0,1),sample_size[1],replace=TRUE,prob=c(py[locus],1-py[locus])))))
        genotype_data <- cbind(genotype_data,geno)
        rm(geno)
      }
      #assign(paste("genotype_data",r,sep="_"),genotype_data)
      
      
      Fstats <- compute_Fstats(genotype_data)
      
      NeHatFc[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFst[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      
      rm(Fstats)
      Fstats <- compute_Fstats(genotype_data,MAF_threshold = 0.1)
      NeHatFcMAF0.1[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFstMAF0.1[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      rm(Fstats)
      Fstats <- compute_Fstats(genotype_data,MAF_threshold = 0.05)
      NeHatFcMAF0.05[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFstMAF0.05[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      rm(Fstats)
      Fstats <- compute_Fstats(genotype_data,MAF_threshold = 0.01)
      NeHatFcMAF0.01[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFstMAF0.01[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      rm(Fstats)
      rm(genotype_data)
    }
    
  }else{# ifelse fis==0
    for ( r in seq_len(number_of_replicates) ){
      if (r==1)     print(paste(simulationID,"Replicate",r))
      if (r%%50==0) print(paste(simulationID,"Replicate",r))
      
      if (model=="Constant") px <- rep(parameter,number_of_loci)   
      if (model=="Beta")     px <- rbeta(number_of_loci,parameter,parameter)    
      id <- paste(model,parameter,Fis,r,sep="_")
      py <- rep(NA,number_of_loci)
      
      
      genotype_data <- cbind( 1:sum(sample_size),
                              c(rep(1,sample_size[1]),rep(2,sample_size[1])) )
      
      for (locus in seq_len(number_of_loci) ){
        f.drift <- px[locus]
        for (g in 1:dT) {
          f.drift <- (rbinom(1,Ne,f.drift)) / Ne
        }
        py[locus] <- f.drift
        
        genotype_prob.x <- c(px[locus]^2 + Fis * (1-px[locus]) * px[locus],
                             2 * (1-px[locus]) * px[locus] * (1-Fis),
                             (1-px[locus])^2 + Fis * (1-px[locus]) * px[locus])
        genotype_prob.y <- c(py[locus]^2 + Fis * (1-py[locus]) * py[locus],
                             2 * (1-py[locus]) * py[locus] * (1-Fis),
                             (1-py[locus])^2 + Fis * (1-py[locus]) * py[locus])
        
        sim_genotypes.x <- as.vector(rmultinom(n    = 1,
                                               size = sample_size[1],
                                               prob = genotype_prob.x))
        sim_genotypes.y <- as.vector(rmultinom(n    = 1,
                                               size = sample_size[1],
                                               prob = genotype_prob.y))
        geno <- c(rep(0,sim_genotypes.x[1]),
                  rep(1,sim_genotypes.x[2]),
                  rep(2,sim_genotypes.x[3]),
                  rep(0,sim_genotypes.y[1]),
                  rep(1,sim_genotypes.y[2]),
                  rep(2,sim_genotypes.y[3]))
        genotype_data <- cbind(genotype_data,geno)
        rm(geno,genotype_prob.x,genotype_prob.y,sim_genotypes.x,sim_genotypes.y)
        
        
        
      }
      #assign(paste("genotype_data",r,sep="_"),genotype_data)
      
      
      Fstats <- compute_Fstats(genotype_data)
      
      NeHatFc[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFst[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      
      rm(Fstats)
      Fstats <- compute_Fstats(genotype_data,MAF_threshold = 0.1)
      NeHatFcMAF0.1[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFstMAF0.1[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      rm(Fstats)
      Fstats <- compute_Fstats(genotype_data,MAF_threshold = 0.05)
      NeHatFcMAF0.05[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFstMAF0.05[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      rm(Fstats)
      Fstats <- compute_Fstats(genotype_data,MAF_threshold = 0.01)
      NeHatFcMAF0.01[r]  <- EstimateNe.F_C(Fstats$F_C,dT,sample_size[1],sample_size[2])
      NeHatFstMAF0.01[r] <- EstimateNe.F_ST(Fstats$F_ST,dT)
      rm(Fstats)
      rm(genotype_data)
    }
  }# ifelse fis==0
  
  
  
  
  boxplot(NeHatFc,NeHatFcMAF0.01,NeHatFcMAF0.05,NeHatFcMAF0.1,
          at=c(1,3,5,7),
          xlim=c(0,9),
          ylim=c(500,1500),
          col="blue",
          #names="Fc",
          outline=F)
  boxplot(NeHatFst,NeHatFstMAF0.01,NeHatFstMAF0.05,NeHatFstMAF0.1,
          at=c(2,4,6,8),
          col="green",
          add=T,
          #names="Fst",
          outline=F)
  abline(h=Ne,col="red")
  
  biasFc         <- mean(NeHatFc         - Ne , na.rm=T)
  biasFst        <- mean(NeHatFst        - Ne , na.rm=T)
  biasFstMAF0.01 <- mean(NeHatFstMAF0.01 - Ne , na.rm=T)
  biasFstMAF0.05 <- mean(NeHatFstMAF0.05 - Ne , na.rm=T)
  biasFstMAF0.1  <- mean(NeHatFstMAF0.1  - Ne , na.rm=T)
  
  maeFc         <- mean(abs(NeHatFc         - Ne) , na.rm=T)
  maeFst        <- mean(abs(NeHatFst        - Ne) , na.rm=T)
  maeFstMAF0.01 <- mean(abs(NeHatFstMAF0.01 - Ne) , na.rm=T)
  maeFstMAF0.05 <- mean(abs(NeHatFstMAF0.05 - Ne) , na.rm=T)
  maeFstMAF0.1  <- mean(abs(NeHatFstMAF0.1  - Ne) , na.rm=T)
  
  mseFc         <- mean((NeHatFc         - Ne)^2 , na.rm=T)
  mseFst        <- mean((NeHatFst        - Ne)^2 , na.rm=T)
  mseFstMAF0.01 <- mean((NeHatFstMAF0.01 - Ne)^2 , na.rm=T)
  mseFstMAF0.05 <- mean((NeHatFstMAF0.05 - Ne)^2 , na.rm=T)
  mseFstMAF0.1  <- mean((NeHatFstMAF0.1  - Ne)^2 , na.rm=T)
  
  rm(f.drift,g,locus,px,py,r)
  save.image(file=paste0("results/SimpleModels/",simulationID,".RData"))
  
} 

