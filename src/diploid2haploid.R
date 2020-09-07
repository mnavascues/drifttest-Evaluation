scenarios <- c(21,23)
num_of_replicates <- 100

load(file="results/simtable.RData")

tau=25

for(scen in scenarios){
  print(paste0("Scenario 0",scen))
  
  Fst_genotypes <- array(NA,num_of_replicates)
  Fst_alleles   <- array(NA,num_of_replicates)

  for(r in 1:num_of_replicates){
    print(paste("Replicate",r))
    
    drifttest_file<- paste0("results/scenario0",scen,"/replicates/scenario0",scen,"_",r,"_drifttest_input")
    diploid_data <- read.table(file=drifttest_file,skip=2)
    
    ncol(diploid_data)
    
    pop1_samples <- which(diploid_data[,1]==1)
    pop2_samples <- which(diploid_data[,1]==2)
    
    n1 <- sum(diploid_data[,1]==1)
    n2 <- sum(diploid_data[,1]==2)
    
    ref_allele_count_pop1 <- as.vector(colSums(diploid_data[pop1_samples,-1]))
    alt_allele_count_pop1 <- n1*2 - ref_allele_count_pop1

    ref_allele_count_pop2 <- as.vector(colSums(diploid_data[pop2_samples,-1]))
    alt_allele_count_pop2 <- n2*2 - ref_allele_count_pop2
    
    haploid_data <- cbind(r1=ref_allele_count_pop1,a1=alt_allele_count_pop1,
                          r2=ref_allele_count_pop2,a2=alt_allele_count_pop2)
    
    tmprldiff_file <- paste0("results/scenario0",scen,"/replicates/scenario0",scen,"_",r,"_tmprldiff_input")
    
    write(2,file=tmprldiff_file,ncolumns=1,append=FALSE)
    write(ncol(diploid_data)-1,file=tmprldiff_file,ncolumns=1,append=TRUE)
    write(t(haploid_data),file=tmprldiff_file,ncolumns=4,append=TRUE)
    
    tmprldiff_commande <- paste("bin/tmprldiff",
                                "-file",tmprldiff_file,
                                "-tau",tau,
                                "-maf 0.000000000001",
                                "-threads 30 > logTmprldiff.txt")
    system(tmprldiff_commande)
    
    res<-read.table(file="logfile.log")
    Fst_alleles[r]<-res[1,1]
    
    file.remove("logfile.log")

    drifttest_commande <- paste("bin/drifttest",
                                "-seed",r,
                                "-tau",tau,
                                "-maf 0.000000000001",
                                "-infile", drifttest_file,
                                "-outfile drifttest_out_",
                                "-threads 30 > logDrifttest.txt")
    system(drifttest_commande)
    
    res<-read.table(file="drifttest_out_multilocus",header=T)
    
    Fst_genotypes[r]<-res[1,1]
    
  }
  assign(paste0("scenario0",scen,"_Fst_genotypes"),Fst_genotypes)
  assign(paste0("scenario0",scen,"_Fst_alleles"),Fst_alleles)
  
    
}
file.remove("outputs.dat")
file.remove("logTmprldiff.txt")
file.remove("logDrifttest.txt")
file.remove("drifttest_out_multilocus")
file.remove("drifttest_out_locus_by_locus")

scenario021_Fst_true <- 0.0123456790123457
scenario023_Fst_true <- 0.0196078431372549

Fst_table<-cbind(out_genot=scenario021_Fst_genotypes,
                 out_allel=scenario021_Fst_alleles,
                 self_genot=scenario023_Fst_genotypes,
                 self_allel=scenario023_Fst_alleles)

boxplot(Fst_table,main="Fst")
lines(x=c(0.5,2.5),y=rep(scenario021_Fst_true,2),col="red",lwd=2)
lines(x=c(2.5,4.5),y=rep(scenario023_Fst_true,2),col="red",lwd=2)

write(t(Fst_table),file="fst_table.txt",ncolumns=4)

