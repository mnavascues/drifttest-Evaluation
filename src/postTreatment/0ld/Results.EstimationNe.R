#######################################################
#
# Evaluation of effective population size estimates
#
#######################################################

source("F_stats_tools.R")

sim_table            <- read.table("Ne/simparams.txt",header=T)
number_of_replicates <- 1000
S1 <- 50
S2 <- 50
total_sample_size    <- 2*(S1+S2)

column_names <- c("id","sigma","deltaT","N","Ne","Fst","NeFst","Fc","NeFc","meanFc","NemeanFc")

results_Ne <- matrix(data = NA,
                     nrow = nrow(sim_table)*number_of_replicates,
                     ncol = length(column_names),
                     dimnames = list(NULL,column_names) )

counter <- 1
for (scenario in seq_along(sim_table$simID) ){
  id     <- sim_table$simID[scenario]
  cat(paste("Scenario",id,"\n"))
  sigma  <- sim_table$sigma[scenario]
  deltaT <- sim_table$selection_period_duration[scenario]
  N      <- sim_table$population_size[scenario]
  Ne     <- (2-sigma)*N
  for (replicate in seq_len(number_of_replicates) ){
    load(paste0("Ne/",id,"/",id,"_",replicate,"_results.RData"))
    Fst      <- results$results_total$F_ST
    NeFst    <- results$results_total$Ne_hat
    Fc       <- results$results_total$F_C
    NeFc     <- EstimateNe.F_C(Fc,deltaT,S1,S2)
    meanFc   <- mean(results$results_by_locus$F_C[results$results_by_locus$maf_test])
    NemeanFc <- EstimateNe.F_C(meanFc,deltaT,S1,S2)
      
    results_Ne[counter,] <- c(scenario,sigma,deltaT,N,Ne,Fst,NeFst,Fc,NeFc,meanFc,NemeanFc)
    counter <- counter+1
  }
  
}




results_Ne <- as.data.frame(results_Ne)

save(results_Ne,sim_table,number_of_replicates,S1,S2,file="Ne/Results.RData")
