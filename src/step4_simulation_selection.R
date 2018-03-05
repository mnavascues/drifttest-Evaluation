# SIMULATION OF THE PERIOD WITH SELECTION (in between samples)

if (!quiet) cat(paste(Sys.time(),"Drift period between replicates for simulation:",simID,", replicate",replic,"\n"))
  
# Drift period between replicates
############################################################################
  
slim_init_drift <- paste0(simID,"_drift_init.txt")
out_drift_lines      <- readLines(con=paste0("bin/",slim_out_drift))
out_drift_pop_line   <- which(out_drift_lines=="Populations:")
write(out_drift_lines[-(1:(out_drift_pop_line-1))], slim_init_drift)

writeMutation          (file=slim_in_drift, number_of_types=2, h=c(0.5,dominance_coef_replic), DFE=c("f","f"), s=c(0,sel_coef_replic), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_drift, u=u)  
writeGenomicElement    (file=slim_in_drift, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
writeChromosome        (file=slim_in_drift, element_type="g1", start=1, end=genome_length)
chromosome_structure <- writeRecombinationChrom(file=slim_in_drift, chr_num=chr_num, genome_length=genome_length, r=r, append=T)
writeGenerations       (file=slim_in_drift, t=drift_period_between_replicates, append=T)
writeDemography        (file=slim_in_drift, type="S", time=1, pop="p1", sigma=sigma)
writeOutput            (file=slim_in_drift, type="A", time=drift_period_between_replicates, filename=slim_out_drift)
writeSeed              (file=slim_in_drift, seed=round(runif(1,-2^31,2^31)) )

write("#INITIALIZATION",file=slim_in_drift,ncolumns=1,append=TRUE)
write(slim_init_drift,file=slim_in_drift,ncolumns=1,append=TRUE)

system(paste("./bin/slim",slim_in_drift,">",slim_log_drift))
system(paste0("mv ",slim_out_drift," bin/",slim_out_drift))
system(paste0("mv ",slim_init_drift," bin/",slim_init_drift))

# An additional generation of drift to add the selective mutation and get first sample
# TO DO: modify here so 1st sample and start of sweep are not simultaneous (additional parameter in sim)
#########################################################################################################

# Set file names
slim_in_drift2  <- paste0("bin/",simID,"_drift2.txt")
slim_out_drift2 <- paste0(simID,"_drift2.out")
slim_log_drift2 <- paste0("bin/",simID,"_drift2.log")

slim_init_drift2 <- paste0(simID,"_drift_init2.txt")
out_drift_lines      <- readLines(con=paste0("bin/",slim_out_drift))
out_drift_pop_line   <- which(out_drift_lines=="Populations:")
write(out_drift_lines[-(1:(out_drift_pop_line-1))], slim_init_drift2)

writeMutation          (file=slim_in_drift2, number_of_types=2, h=c(0.5,dominance_coef_replic), DFE=c("f","f"), s=c(0,sel_coef_replic), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_drift2, u=u)  
writeGenomicElement    (file=slim_in_drift2, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
writeChromosome        (file=slim_in_drift2, element_type="g1", start=1, end=genome_length)
chromosome_structure <- writeRecombinationChrom(file=slim_in_drift2, chr_num=chr_num, genome_length=genome_length, r=r, append=T)
writeGenerations       (file=slim_in_drift2, t=1, append=T)
writeDemography        (file=slim_in_drift2, type="S", time=1, pop="p1", sigma=sigma)
writeOutput            (file=slim_in_drift2, type="A", time=1, filename=slim_out_drift2)
writeSeed              (file=slim_in_drift2, seed=round(runif(1,-2^31,2^31)) )

if (selection_mode=="NM"){
  writePredeterminedMutations(file=slim_in_drift2, time=1, mut_type="m2", x=sample(genome_length,1) )
}

write("#INITIALIZATION",file=slim_in_drift2,ncolumns=1,append=TRUE)
write(slim_init_drift2,file=slim_in_drift2,ncolumns=1,append=TRUE)


system(paste("./bin/slim",slim_in_drift2,">",slim_log_drift2))
system(paste0("mv ",slim_out_drift2," bin/",slim_out_drift2))
system(paste0("mv ",slim_init_drift2," bin/",slim_init_drift2))

# Selection period per se
#################################################

# Set file names
slim_in_selection   <- paste0("bin/",simID,"_selection.txt")
slim_out_selection  <- paste0(simID,"_selection.out")
slim_log_selection  <- paste0("bin/",simID,"_selection.log")
slim_init_selection <- paste0(simID,"_selection_init.txt")

# Read slim output from DRIFT period
out_drift_lines      <- readLines(con=paste0("bin/",slim_out_drift2))
out_drift_pop_line   <- which(out_drift_lines=="Populations:")
out_drift_mut_line   <- which(out_drift_lines=="Mutations:")
out_drift_gen_line   <- which(out_drift_lines=="Genomes:")
out_drift_num_of_mut <- out_drift_gen_line-out_drift_mut_line-1

# Get table of mutations from end of drift period
mutation_table <- read.table(file=paste0("bin/",slim_out_drift2), skip=out_drift_mut_line, nrows=out_drift_num_of_mut)

if (selection_mode!="NM") {
  if (selection_mode!="SV"){
    selection_mode<-"SV"
    warning("Adaptation mode undefined, using standing variation") 
  }
  
  
  
  
  # SELECTION ON STANDING VARIATION (changes selection coefficient of on a random locus)
  
  if (initial_frequency<=0 | initial_frequency>=1){
    # chose a random locus
    sampled_mut <- sample(x=nrow(mutation_table), size=1)
  }else{
    if (advantageous_allele=="derived"){
      initial_allele_count <- round(N*2*initial_frequency)
    }else if (advantageous_allele=="ancestral"){
      initial_allele_count <- round(N*2*(1-initial_frequency))
    }
    if (initial_allele_count==0) initial_allele_count <- 1
    if (initial_allele_count==N*2) initial_allele_count <- N*2-1
    potential_loci <- which(mutation_table[,8]==initial_allele_count)
    lower_limit <- initial_allele_count
    upper_limit <- initial_allele_count
    while (length(potential_loci)==0){
      lower_limit <- lower_limit-1
      if (lower_limit==0) lower_limit <- 1
      upper_limit <- upper_limit+1
      if (upper_limit==N*2) upper_limit <- N*2-1
      potential_loci <- intersect(which(mutation_table[,8]>lower_limit),
                                  which(mutation_table[,8]<upper_limit))
    }
    sampled_mut <- sample(x=potential_loci,size=1)
  }
  
  # change selection coefficient of locus
  levels(mutation_table[,2])    <- c("m1","m2") 
  mutation_table[sampled_mut,2] <- "m2"
  mutation_table[sampled_mut,4] <- sel_coef_replic
  mutation_table[sampled_mut,5] <- dominance_coef_replic
}

# write initialzation file for slim (state of population at starting point of selection period)
write(out_drift_lines[out_drift_pop_line:out_drift_mut_line], slim_init_selection)
write.table(mutation_table, slim_init_selection, append=T, quote=F, col.names=F, row.names=F)
write(out_drift_lines[-(1:(out_drift_gen_line-1))], slim_init_selection, append=T)
#remove(mutation_table)


# write input file for slim
writeMutation          (file=slim_in_selection, number_of_types=2, h=c(0.5,dominance_coef_replic), DFE=c("f","f"), s=c(0,sel_coef_replic), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_selection, u=u)  
writeGenomicElement    (file=slim_in_selection, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
writeChromosome        (file=slim_in_selection, element_type="g1", start=1, end=genome_length)
writeRecombinationChrom(file=slim_in_selection, chr_num=chr_num, genome_length=genome_length, r=r, append=T)  
writeGenerations       (file=slim_in_selection, t=selection_period_duration, append=T)
writeDemography        (file=slim_in_selection, type="S", time=1, pop="p1", sigma=sigma)
writeOutput            (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
writeOutput            (file=slim_in_selection, type="T", time=1, mut_type="m2", append_output=T)
writeOutput            (file=slim_in_selection, type="F", time=selection_period_duration, append_output=T)
writeSeed              (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )

write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)

if (!quiet) cat(paste(Sys.time(),"SIMULATION OF THE PERIOD WITH SELECTION STARTS. Simulation:",simID,"\n"))
system(paste("./bin/slim",slim_in_selection,">",slim_log_selection))
if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF THE PERIOD WITH SELECTION. Simulation:",simID,"\n"))
system(paste0("mv ",slim_out_selection," bin/",slim_out_selection))
system(paste0("mv ",slim_init_selection," bin/",slim_init_selection))
