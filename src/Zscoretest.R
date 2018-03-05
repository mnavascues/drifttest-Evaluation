# run drifttestFSTonly on simulations

project <- "SV"
drifttest_input <- "results/inputD"
num_of_threads <- 10
number_of_scenarios <- 50
number_of_replicates <- 1000
num_of_rep4test <- 100

for (scen in seq_len(number_of_scenarios)){
  
  counter_allele_not_lost <- 0
  replic <- 0
  
  while (counter_allele_not_lost < num_of_rep4test & replic < number_of_replicates) {
    replic <- replic+1  
    
    print(paste("scenario",scen,"; replicate", replic))
    
    RData_results    <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_results.RData")
    RData_data       <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_data.RData")
    RData_params     <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_params.RData")
    RData_trajectory <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_trajectory.RData")
    
    outfile <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_FSTonly_")
    
    load(RData_results)
    
    if (advantageous_allele_not_lost){
      counter_allele_not_lost <- counter_allele_not_lost+1
      if (!file.exists(paste0(outfile,"locus_by_locus"))){
        
        load(RData_data)
        load(RData_params)
        
        number_of_loci <- ncol(SNP_data$genotype_data)-2
        number_of_ind  <- nrow(SNP_data$genotype_data)
        
        write(number_of_ind,file=drifttest_input,ncolumns=1,append=F)
        write(number_of_loci,file=drifttest_input,ncolumns=1,append=T)
        
        write.table(SNP_data$genotype_data[,-1],
                    file=drifttest_input,
                    row.names=F,
                    col.names=F,
                    append=T)
        
        system(paste("./bin/drifttestFSTonly",
                     "-tau", selection_period_duration,
                     "-seed", round(runif(1,0,10000)),
                     "-threads",num_of_threads,
                     "-maf",MAF_threshold,
                     "-infile", drifttest_input,
                     "-outfile",outfile) )
        
        file.remove(drifttest_input)
        
      }else{
        print(paste("Drifttest already run for replicate",replic,": skipping drifttest"))
      }
      
    }else{
      print(paste("Advantageous allele lost in replicate",replic,": skipping drifttest"))
    }
    
    
    
  }
}

