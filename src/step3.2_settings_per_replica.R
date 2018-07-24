# file names for replicates

whole_pop_file    <- paste0(working_dir,"/replicates/",simID,"_",replic,"_whole_pop.RData") 
data_file         <- paste0(working_dir,"/replicates/",simID,"_",replic,"_data.RData") 
drifttest_infile  <- paste0(working_dir,"/replicates/",simID,"_",replic,"_drifttest_input") 
drifttest_outfile <- paste0(working_dir,"/replicates/",simID,"_",replic,"_drifttest_out_") 

if (selection_mode=="NM"){
  sel_coef_replic       <- sel_coef
  dominance_coef_replic <- dominance_coef
  
}else{
  if (selection_mode!="SV"){
    selection_mode<-"SV"
    warning("Adaptation mode undefined, using standing variation") 
  }
  # choose advantageous allele between derived and ancestral state
  if(is.na(fixed_advantageous_allele)){
    advantageous_allele <- sample(c("derived","ancestral"),size=1)
  }else if (fixed_advantageous_allele=="derived" | fixed_advantageous_allele=="ancestral"){
    advantageous_allele <- fixed_advantageous_allele
  }else{
    advantageous_allele <- sample(c("derived","ancestral"),size=1)
  }
  
  if (advantageous_allele=="derived"){
    sel_coef_replic       <- sel_coef
    dominance_coef_replic <- dominance_coef
    write(paste("Derived allele is advantageous in replicate",replic), file=log_file,append=T)
  }else if (advantageous_allele=="ancestral"){
    sel_coef_replic       <- -(sel_coef/(1+sel_coef))
    dominance_coef_replic <- 1 - dominance_coef
    write(paste("Ancestral allele is advantageous in replicate",replic), file=log_file,append=T)
  }
}

