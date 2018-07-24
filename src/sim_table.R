number_of_sims <- 60

sim_table <- data.frame(simID=array(NA,number_of_sims),
                        sigma=array(NA,number_of_sims),
                        u=array(NA,number_of_sims),
                        r=array(NA,number_of_sims),
                        genome_length=array(NA,number_of_sims),
                        N=array(NA,number_of_sims),
                        sel_coef=array(NA,number_of_sims),
                        dominance_coef=array(NA,number_of_sims),
                        selection_mode=array(NA,number_of_sims),
                        initial_frequency=array(NA,number_of_sims),
                        selection_period_duration=array(NA,number_of_sims),
                        sample_size_t1=array(NA,number_of_sims),
                        sample_size_t2=array(NA,number_of_sims),
                        sample_size_loci=array(NA,number_of_sims),
                        MAF_threshold=array(NA,number_of_sims),
                        number_of_replicates=array(NA,number_of_sims))

for (sim in 1:number_of_sims){
  
  if (sim<10)       simID <- paste0("scenario00",sim)
  else if (sim<100) simID <- paste0("scenario0",sim)
  else              simID <- paste0("scenario",sim)
  
  load(paste0("results/",simID,"/",simID,"_params.RData"))
  
  sim_table$simID[sim]                     <- simID
  sim_table$sigma[sim]                     <- sigma
  sim_table$u[sim]                         <- u
  sim_table$r[sim]                         <- r
  sim_table$genome_length[sim]             <- genome_length
  sim_table$N[sim]                         <- N
  sim_table$sel_coef[sim]                  <- sel_coef
  sim_table$dominance_coef[sim]            <- dominance_coef
  sim_table$selection_mode[sim]            <- selection_mode
  sim_table$initial_frequency[sim]         <- initial_frequency
  sim_table$selection_period_duration[sim] <- selection_period_duration
  sim_table$sample_size_t1[sim]            <- sample_size[1]
  sim_table$sample_size_t2[sim]            <- sample_size[2]
  sim_table$sample_size_loci[sim]          <- sample_size_loci
  sim_table$MAF_threshold[sim]             <- MAF_threshold
  sim_table$number_of_replicates[sim]      <- number_of_replicates

}
save(sim_table,file="results/simtable.RData")


