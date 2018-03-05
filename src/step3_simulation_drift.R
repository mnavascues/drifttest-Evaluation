# SIMULATION OF A PURE DRIFT PERIOD


# Set file names
slim_in_drift  <- paste0("bin/",simID,"_drift.txt")
slim_out_drift <- paste0(simID,"_drift.out")
slim_log_drift <- paste0("bin/",simID,"_drift.log")

# Write slim input file (functions from slim_tools_2.R)
writeMutation          (file=slim_in_drift, number_of_types=2, h=c(0.5,dominance_coef), DFE=c("f","f"), s=c(0,sel_coef), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_drift, u=u)  
writeGenomicElement    (file=slim_in_drift, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
writeChromosome        (file=slim_in_drift, element_type="g1", start=1, end=genome_length)
chromosome_structure <- writeRecombinationChrom(file=slim_in_drift, chr_num=chr_num, genome_length=genome_length, r=r, append=T)
writeGenerations       (file=slim_in_drift, t=drift_period_duration, append=T)
if (file.exists(paste0("results/",simID,"/populationATequilibrium"))){
  writeDemography        (file=slim_in_drift, type="S", time=1, pop="p1", sigma=sigma)
}else{
  writeDemography        (file=slim_in_drift, type="P", time=1, pop="p1", N=N)
  writeDemography        (file=slim_in_drift, type="S", time=1, pop="p1", sigma=sigma, append_demography=T)
}
writeOutput            (file=slim_in_drift, type="A", time=drift_period_duration, filename=slim_out_drift)
writeSeed              (file=slim_in_drift, seed=round(runif(1,-2^31,2^31)) )
if (file.exists(paste0("results/",simID,"/populationATequilibrium"))){
  file.copy(from = paste0("results/",simID,"/populationATequilibrium"),
            to   = "populationATequilibrium",
            overwrite = TRUE)
  write("#INITIALIZATION",file=slim_in_drift,ncolumns=1,append=TRUE)
  write("populationATequilibrium",file=slim_in_drift,ncolumns=1,append=TRUE)
}

if (!quiet) cat(paste(Sys.time(),"SIMULATION OF A PURE DRIFT PERIOD STARTS. Simulation:",simID, "\n"))

# Run SLiM
system(paste("./bin/slim",slim_in_drift,">",slim_log_drift))
# NB: SLiM executable must be in the same folder as SLiM input file (this is a requirement from SLiM)
system(paste0("mv ",slim_out_drift," bin/",slim_out_drift))

if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF A PURE DRIFT PERIOD. Simulation:",simID,"\n"))


