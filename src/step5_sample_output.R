
# Read slim output and log from SELECTION period
#------------------------------------------------

out_selection_lines    <- readLines(con=paste0("bin/",slim_out_selection))
out_selection_mut_line <- which(out_selection_lines=="Mutations:")
out_selection_gen_line <- which(out_selection_lines=="Genomes:")

log_selection_lines    <- readLines(con=slim_log_selection)
log_selection_mut_line <- which(log_selection_lines=="Mutations:")
if (length(log_selection_mut_line)>0){
  log_selection_num_mut <- length(grep(pattern = "#OUT",
                                       x = log_selection_lines[seq_along(log_selection_lines)>log_selection_mut_line],
                                       invert = TRUE))
}else{
  log_selection_num_mut <- 0
}

# VERIFYING THE SUCCESS OF SELECTION (was the advantageous allele lost by drift?)
#---------------------------------------------------------------------------------

if (advantageous_allele=="derived"){
  if (log_selection_num_mut > 0) {
    if (length(grep(pattern="m2",x=out_selection_lines[(out_selection_mut_line+1):(out_selection_gen_line-1)]))>0 || length(grep(pattern="m2",x=log_selection_lines[(log_selection_mut_line+1):length(log_selection_lines)]))>0){
      advantageous_allele_not_lost <- TRUE
    }else{ advantageous_allele_not_lost <- FALSE }  
  }else if(log_selection_num_mut==0){
    if (length(grep(pattern="m2",x=out_selection_lines[(out_selection_mut_line+1):(out_selection_gen_line-1)]))>0){
      advantageous_allele_not_lost <- TRUE
    }else{ advantageous_allele_not_lost <- FALSE }
  }
}else if (advantageous_allele=="ancestral"){
  if (log_selection_num_mut > 0) {
    if (length(grep(pattern="m2",x=log_selection_lines[log_selection_mut_line+1:log_selection_num_mut]))>0){
      advantageous_allele_not_lost <- FALSE
    }else{ advantageous_allele_not_lost <- TRUE }  
  }else if(log_selection_num_mut==0){
    advantageous_allele_not_lost <- TRUE
  }
}
if (skip_selection) advantageous_allele_not_lost <- FALSE

if(!advantageous_allele_not_lost){
  if (!quiet) cat(paste(Sys.time(),"Advantageous allele was lost in simulation",simID,"replicate",replic, "\n"))
  #write(paste("Advantageous allele was lost in replicate",replic), file=log_file,append=T)
}else{
  if (!quiet) cat(paste(Sys.time(),"Advantageous allele was NOT lost in simulation",simID,"replicate",replic, "\n"))
  #write(paste("Advantageous allele was NOT lost in replicate",replic), file=log_file,append=T)
}  


if (advantageous_allele_not_lost | sel_coef==0){
  
  # SAVING LOCUS UNDER SELECTION ALLELE FREQUENCY TRAJECTORY
  #----------------------------------------------------------
  trajectory <- log_selection_lines[grep(pattern = "#OUT:",x = log_selection_lines)]
  trajectory <- trajectory[grep(pattern = " T ",x = trajectory)]
  if (length(trajectory)>0){
    trajectory <- as.numeric(matrix(unlist(strsplit(trajectory," ")),ncol =11,byrow=T)[,11])
  }
  
  
  # Make a data frame with ALL polymorphic sites, and matrix of haplotypes
  #-----------------------------------------------------------------------
  SNP_table <- Make_SNP_table(file_t1     = paste0("bin/",slim_init_selection),
                              file_t2     = paste0("bin/",slim_out_selection),
                              file_fixed  = slim_log_selection)
  
  if (!quiet) cat(paste(Sys.time(),"List of all polymorphic SNP completed for simulation",simID,"replicate",replic, "\n"))     
  haplotypes_1 <- SNP_table$haplotypes_1
  haplotypes_2 <- SNP_table$haplotypes_2
  m2           <- SNP_table$m2
  removed_loci <- unique(SNP_table$removed_loci)
  SNP_table    <- SNP_table$SNP_table
  
  
  total_S <- nrow(SNP_table)
  write(paste("Total number of bi-allelic loci in the replicate",replic,":",total_S), file=log_file,append=T)

  # Make a sample of the number of individuals and number of loci specified
  #-------------------------------------------------------------------------
  SNP_data <- Make_sample(haplotypes_1,
                          haplotypes_2,
                          SNP_table,
                          removed_loci,
                          sample_size,
                          sample_size_loci,
                          Ne_only)
  
  if (!no_whole_pop_out){
    whole_pop_data <-  Make_sample(haplotypes_1,
                                   haplotypes_2,
                                   SNP_table,
                                   removed_loci,
                                   c(N,N),
                                   total_S,
                                   Ne_only)
    save(whole_pop_data,m2,SNP_table,file=whole_pop_file)
    write(paste("Whole population configuration (SNP table + haplotypes) saved into:",whole_pop_file), file=log_file,append=T)
  }
  
  
  if (!Ne_only){
  # write drifttest input file  
  write(sum(sample_size), file=drifttest_infile)
  write(sample_size_loci, file=drifttest_infile,append=T)
  write(t(SNP_data$genotype_data[,-1]),
        file=drifttest_infile,
        ncolumns = sample_size_loci+1 ,
        append=T)
  }
}



