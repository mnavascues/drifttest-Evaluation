######################################################
#
# Run neutrality test on already existing simulations
#
######################################################
library(qqman)
source("F_stats_tools.R")

sim_table            <- read.table("Ne/simparams.txt",header=T)
number_of_replicates <- 1

MAF_threshold   <- 0.01
num_of_sim_test <- 100
Ne_only         <- FALSE

# choose set of simulations
selection_period_duration <- 25
population_size           <- 500

scenarios <- intersect(which(sim_table$selection_period_duration==selection_period_duration),
                       which(sim_table$population_size==population_size))

for (scen in 10:10){
  print(paste("Scenario",scenarios[scen]))
  for (rep in 4:5){
    print(paste("Replicate",rep))
    

    load(file=paste0("Ne/Ne",scenarios[scen],"/Ne",scenarios[scen],"_",rep,"_data.RData"))
    
    
    pdf(file=paste0("SLiM_Ne",scenarios[scen],"_",rep,"_QQplot.pdf"),paper="a4")
    # DIRICHLET
    results_file <- paste0("Ne/Ne",scenarios[scen],"/Ne",scenarios[scen],"_",rep,"_resultsDirichlet.RData") 
    results <- FST_outlier_test.Fis(SNP_data$genotype_data,
                                    MAF_threshold=MAF_threshold,
                                    delta_T=selection_period_duration,
                                    num_of_sim_test=num_of_sim_test)
    qqman::qq(results$results_by_locus$p_value,type="l",xlim=c(0,5),ylim=c(0,5))
    save(results,file=results_file)
    # BETA
    results_file <- paste0("Ne/Ne",scenarios[scen],"/Ne",scenarios[scen],"_",rep,"_resultsBeta.RData") 
    results <- FST_outlier_test(SNP_data$genotype_data,
                                    MAF_threshold=MAF_threshold,
                                    delta_T=selection_period_duration,
                                    num_of_sim_test=num_of_sim_test)
    qqman::qq(results$results_by_locus$p_value,type="l",xlim=c(0,5),ylim=c(0,5))
    save(results,file=results_file)
    dev.off()
    

    #chromosome                                      <- array(1,length(SNP_data$SNP_table$x))
    #chromosome[which(SNP_data$SNP_table$x>(5e8/2))] <- 2

    #position                <- SNP_data$SNP_table$x
    #position[chromosome==2] <- position[chromosome==2]-(5e8/2)
    
    #res <- data.frame(SNP=SNP_data$SNP_table$type,CHR=chromosome,BP=position,P=results$results_by_locus$p_value)
    #res <- res[results$results_by_locus$maf_test,]
    #if(!is.na(results$results_total$Ne_hat)){
    #  manhattan(res,ylim=c(0,5))
    #  qqman::qq(res$P,type="l",xlim=c(0,5),ylim=c(0,5))
    #}
  }
}














kk <- data.frame(ks.D1    = array(NA,length(scenarios)*number_of_replicates),
                 ks.D2    = array(NA,length(scenarios)*number_of_replicates),
                 ks.D3    = array(NA,length(scenarios)*number_of_replicates),
                 Ne_error = array(NA,length(scenarios)*number_of_replicates),
                 Ne_bias = array(NA,length(scenarios)*number_of_replicates),
                 Ne_hat   = array(NA,length(scenarios)*number_of_replicates),
                 sigma    = array(NA,length(scenarios)*number_of_replicates))

counter <- 1
for (scen in seq_along(scenarios)){
  print(paste("Scenario",scenarios[scen]))
  
  for (rep in 1:number_of_replicates){
    print(paste("Replicate",rep))
    load(file=paste0("Ne/Ne",scenarios[scen],"/Ne",scenarios[scen],"_",rep,"_results2.RData"))

    sigma <- sim_table[sim_table$simID==paste0("Ne",scenarios[scen]),"sigma"]
    N     <- sim_table[sim_table$simID==paste0("Ne",scenarios[scen]),"population_size"]
    Ne    <- (2-sigma)*N
    
    kk$ks.D1[counter]     <- ks.test(results$results_by_locus$p_value,punif)$statistic
    kk$Ne_hat[counter]   <- results$results_total$Ne_hat
    kk$Ne_error[counter] <- sqrt((results$results_total$Ne_hat-Ne)^2)/Ne
    kk$Ne_bias[counter]  <- (results$results_total$Ne_hat-Ne)/Ne
    kk$sigma[counter]    <- sigma

    load(file=paste0("Ne/Ne",scenarios[scen],"/Ne",scenarios[scen],"_",rep,"_results3.RData"))
    kk$ks.D2[counter]     <- ks.test(results$results_by_locus$p_value,punif)$statistic
    load(file=paste0("Ne/Ne",scenarios[scen],"/Ne",scenarios[scen],"_",rep,"_results4.RData"))
    kk$ks.D3[counter]     <- ks.test(results$results_by_locus$p_value,punif)$statistic
    
    counter <- counter+1
  }
}




plot(kk$Ne_hat,kk$ks.D1)
plot(kk$Ne_hat,kk$ks.D2)
#plot(kk$Ne_hat,kk$ks.D3)

plot(kk$Ne_bias,kk$ks.D1)
plot(kk$Ne_bias,kk$ks.D2)
plot(kk$Ne_bias,kk$ks.D3)

plot(kk$Ne_bias,kk$ks.D1)
for (scen in seq_along(scenarios)){
  
  i <- 1:100+(scen-1)*100

  plot(kk$Ne_bias[i],kk$ks.D1[i], main=scenarios[scen],ylim=c(0,0.5))
  reg<-lm(kk$ks.D1[i]~kk$Ne_bias[i])
  abline(reg, col="red")
  
  #plot(kk$Ne_bias[i],kk$ks.D2[i],main=scenarios[scen])
  #reg<-lm(kk$ks.D2[i]~kk$Ne_bias[i])
  #abline(reg, col="red")

  plot(kk$Ne_bias[i],kk$ks.D3[i],main=scenarios[scen],ylim=c(0,0.5))
  reg<-lm(kk$ks.D3[i]~kk$Ne_bias[i])
  abline(reg, col="red")
  
}
#plot(kk$Ne_bias,kk$ks.D3)

plot(kk$sigma,kk$ks.D1)
boxplot(kk$ks.D1,kk$ks.D3)

boxplot(kk$ks.D1~kk$sigma,ylim=c(0,0.5))
boxplot(kk$ks.D3~kk$sigma,ylim=c(0,0.5))
