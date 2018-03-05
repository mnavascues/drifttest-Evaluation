require(gtools)

#	The input data file should be formatted as follows:
#
#	1	1	0	1	2	2	1	0
#	2	1	0	1	0	0	1	2
#	1	2	1	1	0	0	1	2
#	2	2	2	1	1	2	0	2
#	3	2	2	2	1	1	1	1
#
#	Each line corresponds to an individual
# 	The first column contains IDs for individual
# 	The second column contains IDs for sampled demes
# 	The next columns correspond to loci
# 	`0' corresponds to homozygotes for type-1 alleles
# 	`1' corresponds to heterozygotes
# 	`2' corresponds to homozygotes for type-2 alleles

# data as an R matrix with appropriate format or infile as a text file in the appropriate format
# MAF_threshold:   threshold to filter loci with a minumim maf (minor allele frequency)
# delta_T:         number of generations between time samples
# num_of_sim_test: number of replicates of drift simulations to build null FC distribution
FST_outlier_test <- function(data,infile=NA,MAF_threshold,delta_T,num_of_sim_test,Fst=NA){
  # read the data
  if (!is.na(infile)) data <- read.table(infile) 
  
  nbr.pops <- length(unique(data[,2]))																# compute the number of samples
  nbr.loci <- ncol(data)-2																						# compute the number of loci
  stopifnot(nbr.pops==2) 
  
  # compute Fst, Fis, Fc
  Fstats   <- compute_Fstats(data,MAF_threshold=MAF_threshold)
  Ne_hat_FST  <- EstimateNe.F_ST(Fstats$F_ST,delta_T)
  Ne_hat      <- round(Ne_hat_FST)

  p_value <- array(NA,dim=nbr.loci)
  if(is.na(Ne_hat)){
    results_total    <- list(F_ST  =Fstats$F_ST,
                             F_IS  =Fstats$F_IS,
                             F_C   =Fstats$F_C,
                             Ne_hat=Ne_hat_FST)
    results_by_locus <- data.frame(F_ST=Fstats$F_ST_locus,
                                   F_IS=Fstats$F_IS_locus,
                                   F_C =Fstats$F_C_locus,
                                   p_value,
                                   maf=Fstats$maf$maf,
                                   maf_test=Fstats$maf$test)
    
  }else{
    for (locus in seq_len(nbr.loci)) {
      if(Fstats$maf$test[locus]){
        genotypes_locus <- rbind(c(sum(data[data[,2]==1,2+locus]==0),
                                   sum(data[data[,2]==1,2+locus]==1),
                                   sum(data[data[,2]==1,2+locus]==2)),
                                 c(sum(data[data[,2]==2,2+locus]==0),
                                   sum(data[data[,2]==2,2+locus]==1),
                                   sum(data[data[,2]==2,2+locus]==2)))
        p_value[locus] <- Drift.simulation.4.test(genotypes = genotypes_locus,
                                                  Ne  = Ne_hat,
                                                  dT  = delta_T,
                                                  num_of_sim_test = num_of_sim_test,
                                                  MAF_threshold   = MAF_threshold)
      }
    }
    
    results_total    <- list(F_ST  =Fstats$F_ST,
                             F_IS  =Fstats$F_IS,
                             F_C   =Fstats$F_C,
                             Ne_hat=Ne_hat_FST)
    if (nbr.loci==1){
      results_by_locus <- data.frame(F_ST=Fstats$F_ST,
                                     F_IS=Fstats$F_IS,
                                     F_C =Fstats$F_C,
                                     p_value,
                                     maf=Fstats$maf$maf,
                                     maf_test=Fstats$maf$test)
    }else{
      results_by_locus <- data.frame(F_ST=Fstats$F_ST_locus,
                                     F_IS=Fstats$F_IS_locus,
                                     F_C =Fstats$F_C_locus,
                                     p_value,
                                     maf=Fstats$maf$maf,
                                     maf_test=Fstats$maf$test)
    }
  }
  return(list(results_total=results_total,results_by_locus=results_by_locus)) 
}


FST_outlier_test.Fis <- function(data,infile=NA,MAF_threshold,delta_T,num_of_sim_test){
  # read the data
  if (!is.na(infile)) data <- read.table(infile) 
  
  nbr.pops <- length(unique(data[,2]))																# compute the number of samples
  nbr.loci <- ncol(data)-2																						# compute the number of loci
  stopifnot(nbr.pops==2) 
  
  # compute Fst, Fis, Fc
  Fstats     <- compute_Fstats(data,MAF_threshold=MAF_threshold)
  Ne_hat_FST <- EstimateNe.F_ST(Fstats$F_ST,delta_T)
  Ne_hat     <- round(Ne_hat_FST)

  p_value <- array(NA,dim=nbr.loci)
  if(is.na(Ne_hat)){
    results_total    <- list(F_ST  =Fstats$F_ST,
                             F_IS  =Fstats$F_IS,
                             F_C   =Fstats$F_C,
                             Ne_hat=Ne_hat_FST)
    results_by_locus <- data.frame(F_ST=Fstats$F_ST_locus,
                                   F_IS=Fstats$F_IS_locus,
                                   F_C =Fstats$F_C_locus,
                                   p_value,
                                   maf=Fstats$maf$maf,
                                   maf_test=Fstats$maf$test)
  }else{
    for (locus in seq_len(nbr.loci)) {
      if(Fstats$maf$test[locus]){
        genotypes_locus <- rbind(c(sum(data[data[,2]==1,2+locus]==0),
                                   sum(data[data[,2]==1,2+locus]==1),
                                   sum(data[data[,2]==1,2+locus]==2)),
                                 c(sum(data[data[,2]==2,2+locus]==0),
                                   sum(data[data[,2]==2,2+locus]==1),
                                   sum(data[data[,2]==2,2+locus]==2)))
        p_value[locus] <- Drift.simulation.4.test.Fis(genotypes = genotypes_locus,
                                                      Ne  = Ne_hat,
                                                      dT  = delta_T,
                                                      Fis = Fstats$F_IS,
                                                      num_of_sim_test = num_of_sim_test,
                                                      MAF_threshold   = MAF_threshold)
      }
    }
    
    results_total    <- list(F_ST  =Fstats$F_ST,
                             F_IS  =Fstats$F_IS,
                             F_C   =Fstats$F_C,
                             Ne_hat=Ne_hat_FST)
    if (nbr.loci==1){
      results_by_locus <- data.frame(F_ST=Fstats$F_ST,
                                     F_IS=Fstats$F_IS,
                                     F_C =Fstats$F_C,
                                     p_value,
                                     maf=Fstats$maf$maf,
                                     maf_test=Fstats$maf$test)
    }else{
      results_by_locus <- data.frame(F_ST=Fstats$F_ST_locus,
                                     F_IS=Fstats$F_IS_locus,
                                     F_C =Fstats$F_C_locus,
                                     p_value,
                                     maf=Fstats$maf$maf,
                                     maf_test=Fstats$maf$test)
    }
  }
  return(list(results_total=results_total,results_by_locus=results_by_locus)) 
}

# Simulates drift from an initial allele frequency, an effective population size and
# a given number of generations to obtain the neutral distribution of the Fst statistic
Drift.simulation.4.test.Fis <- function(genotypes,Ne,dT,Fis,num_of_sim_test,MAF_threshold){
  
  # require(gtools)
  # genotypes<-rbind(c(23,5,22),c(18,2,30))
  
  if(Fis<0) Fis <- 0
  Ne  <- round(Ne)
  Fst <- compute_FST_from_genotype_counts(genotypes)
  sample_size <- rowSums(genotypes)
  pvalue <- 0
  sim <- 1
  while (sim <= num_of_sim_test) {
    genotype_prob.x <- rdirichlet(n     = 1,
                                  alpha = c(1+genotypes[1,1],1+genotypes[1,2],1+genotypes[1,3]))
    sim_genotypes.x <- as.vector(rmultinom(n    = 1,
                                           size = sample_size[1],
                                           prob = genotype_prob.x))
    x  <- sum( sim_genotypes.x * c(2,1,0) )
    px <- genotype_prob.x[1]+genotype_prob.x[2]/2
    py <- px
    for (g in 1:dT) {
      py <- (rbinom(1,Ne,py)) / Ne
    }
    genotype_prob.y <- c(py^2 + Fis * (1-py) * py,
                         2 * (1-py) * py * (1-Fis),
                         (1-py)^2 + Fis * (1-py) * py)
    sim_genotypes.y <- as.vector(rmultinom(n    = 1,
                                           size = sample_size[2],
                                           prob = genotype_prob.y))
    y  <- sum( sim_genotypes.y * c(2,1,0) )
    
    mean.freq <- mean(c(x,y)/(sample_size*2))
    
    sim_genotypes <- rbind(sim_genotypes.x,sim_genotypes.y)
    if ( mean.freq >= MAF_threshold & mean.freq <= (1-MAF_threshold) ) {
      Fst_sim <- compute_FST_from_genotype_counts(sim_genotypes)
      if (Fst_sim >= Fst) pvalue <- pvalue + 1
      sim <- sim + 1
    }
    if ((sim == num_of_sim_test) & (pvalue == 0) & (num_of_sim_test <= 1e4)) {
      num_of_sim_test <- num_of_sim_test*10
    }
  }
  pvalue <- pvalue/num_of_sim_test
  return(pvalue)
}


# Simulates drift from an initial allele frequency, an effective population size and
# a given number of generations to obtain the neutral distribution of the Fc statistic
Drift.simulation.4.test <- function(genotypes,Ne,dT,num_of_sim_test,MAF_threshold){
  
  # require(gtools)
  # genotypes<-rbind(c(23,5,22),c(18,2,30))
  
  Ne          <- round(Ne)
  Fst         <- compute_FST_from_genotype_counts(genotypes)
  sample_size <- rowSums(genotypes)
  n           <- sample_size[1]*2 
  k           <- genotypes[1,1]*2+genotypes[1,2]

  pvalue <- 0
  sim <- 1
  while (sim <= num_of_sim_test) {
    px          <- rbeta( n=1, shape1=1+k, shape2=1+n-k )
    x           <- rbinom(n=1, size=n,prob=px)
    py <- px
    for (g in 1:dT) {
      py <- (rbinom(1,Ne,py)) / Ne
    }
    y <- rbinom(1, 2*sample_size[2],py)
    mean.freq <- mean(c(x,y)/(sample_size*2))
    if ( mean.freq >= MAF_threshold & mean.freq <= (1-MAF_threshold) ) {
      genotypes.x <- sample(c(rep(0,x),rep(1,n-x)))
      genotypes.x <- rbind(genotypes.x[1:(length(genotypes.x)/2)],
                           genotypes.x[(length(genotypes.x)/2+1):length(genotypes.x)])
      genotypes.x <- colSums(genotypes.x)
      genotypes.x <- c(sum(genotypes.x==0),sum(genotypes.x==1),sum(genotypes.x==2))
      genotypes.y<-sample(c(rep(0,y),rep(1,2*sample_size[2]-y)))
      genotypes.y <- rbind(genotypes.y[1:(length(genotypes.y)/2)],
                           genotypes.y[(length(genotypes.y)/2+1):length(genotypes.y)])
      genotypes.y <- colSums(genotypes.y)
      genotypes.y <- c(sum(genotypes.y==0),sum(genotypes.y==1),sum(genotypes.y==2))
      sim_genotypes <- rbind(genotypes.x,genotypes.y)

      Fst_sim <- compute_FST_from_genotype_counts(sim_genotypes)
      if (Fst_sim >= Fst) pvalue <- pvalue + 1
      sim <- sim + 1
    }
    if ((sim == num_of_sim_test) & (pvalue == 0) & (num_of_sim_test <= 1e4)) {
      num_of_sim_test <- num_of_sim_test*10
    }
  }
  pvalue <- pvalue/num_of_sim_test
  return(pvalue)
}





# Calculates Kimura probability density
dens_kim<-function(f,p,c,eps=1e-8){ #ici c=t/2N avec 2N= nombre efficace haploides
  c=c*2 #dans le programme c'est t/N (on multiplie apres pra 0.25)
  f=min(1,max(0,f))
  ##cas fixe à 0
  if(f==0){
    r=1-2*(1-p) ; niter=2 ; Tr_0=1 ; Tr_1=3*r
    dens= -1.5*Tr_0*exp(-0.5*c) + (5/6)*Tr_1*exp(-1.5*c)
    crit=0 ; stop=0
    while(crit!=2){ #pednant dix iterations on tmpd negligeable
      niter=niter+1 ; tmpn=niter-1 #on regarde a i-1 les polynomes
      tmpr=(2*r*(tmpn+0.5)*Tr_1 - (tmpn+1)*Tr_0)/tmpn
      tmpi=(2*niter+1)/(niter*(niter+1)) ; tmpd=(exp(log(tmpi)-0.25 *niter*(niter+1)*c))*tmpr
      if(abs(tmpd)<eps){crit=crit+1}else{crit=0}
      if(niter%%2==0){dens=dens+tmpd}else{dens=dens-tmpd}
      Tr_0=Tr_1 ; Tr_1=tmpr
    }
    dens=(1-r**2)*dens/2 + 1 - p
  }
  #cas fixe à 1
  if(f==1){
    r=1-2*p ; niter=2 ; Tr_0=1 ; Tr_1=3*r
    dens=-1.5*Tr_0*exp(-0.5*c) + (5/6)*Tr_1*exp(-1.5*c)
    crit=0 ; stop=0
    while(crit!=2){ #pednant dix iterations on tmpd negligeable
      niter=niter+1 ; tmpn=niter-1 #on regarde a i-1 les polynomes
      tmpr=(2*r*(tmpn+0.5)*Tr_1 - (tmpn+1)*Tr_0)/tmpn
      tmpi=(2*niter+1)/(niter*(niter+1)) ; tmpd=(exp(log(tmpi)-0.25 *niter*(niter+1)*c))*tmpr
      if(abs(tmpd)<eps){crit=crit+1}else{crit=0}
      if(niter%%2==0){dens=dens+tmpd}else{dens=dens-tmpd}
      Tr_0=Tr_1  ; Tr_1=tmpr
    }
    dens=(1-r**2)*dens/2 + p
  }
  ##cas non fixe
  if(f>0 & f<1){
    r=1-2*p ; z=1-2*f
    niter=2
    Tr_0=Tz_0=1 #on stocke T(n-2) la dedans
    Tr_1=3*r ; Tz_1=3*z  #on strocke T(n-1) la dedans
    dens=1.5*Tr_0*Tz_0*exp(-0.5*c) + (5/6)*Tr_1*Tz_1*exp(-1.5*c)
    crit=0 ; stop=0
    while(crit!=2){ #pednant dix iterations on tmpd negligeable
      niter=niter+1 #terme de la somme
      tmpn=niter-1 #on regarde a i-1 les polynomes
      tmpr=(2*r*(tmpn+0.5)*Tr_1 - (tmpn+1)*Tr_0)/tmpn  ; tmpz=(2*z*(tmpn+0.5)*Tz_1 - (tmpn+1)*Tz_0)/tmpn
      tmpi=(2*niter+1)/(niter*(niter+1)) ; tmpd=(exp(log(tmpi)-0.25 *niter*(niter+1)*c))*tmpr*tmpz
      if(abs(tmpd)<eps){crit=crit+1}else{crit=0}
      dens=dens+tmpd
      Tr_0=Tr_1 ; Tz_0=Tz_1 ; Tr_1=tmpr ; Tz_1=tmpz
    }
    dens=(1-r**2)*dens}
  
  list(dens=dens,niter=niter)
}


# Estimates p-value from Kimura probability density
Kimura.density.4.test <- function(Ne,dT,starting_freq,sample_size,grid_size,MAF_threshold){
}



####### Create filter for loci using a maf threshold
maf_filter <- function(p,MAF_threshold){
  maf          <- (p[1,]+p[2,])/2
  maf[maf>0.5] <- 1- maf[maf>0.5]
  test         <- maf >= MAF_threshold
  return(list(maf=maf,test=test))
}





######## Estimating Effective Population Size (Ne)
EstimateNe.F_C <- function(FC,dT,S1,S2) {
  Ne_hat <- 2*((dT -2)/ (2*(FC - (1/(2*S1)) - (1/(2*S2)))))
  if (Ne_hat<0) Ne_hat <- NA
  return(Ne_hat)
}
EstimateNe.F_ST <- function(Fst_hat,dT) {
  Ne_hat <- dT * (1 - Fst_hat) / (2 * Fst_hat) 
  if (Ne_hat<0) Ne_hat <- NA
  return(Ne_hat)
}


#	The input data file should be formatted as follows:
#
#	1	1	0	1	2	2	1	0
#	2	1	0	1	0	0	1	2
#	1	2	1	1	0	0	1	2
#	2	2	2	1	1	2	0	2
#	3	2	2	2	1	1	1	1
#
#	Each line corresponds to an individual
# 	The first column contains IDs for individual
# 	The second column contains IDs for sampled demes
# 	The next columns correspond to loci
# 	`0' corresponds to homozygotes for type-1 alleles
# 	`1' corresponds to heterozygotes
# 	`2' corresponds to homozygotes for type-2 alleles
#
#	In the above example, there are two sampled demes made of 2 and 3 diploid individuals, respectively
#	For the above example, Genepop gives the following estimates of F-statistics:
#
#	Multilocus estimates for diploid data
#	Locus           Fwc(is)     Fwc(st)     Fwc(it)
#	------------    -------     -------     -------
#	loc_1            0.0526      0.7543      0.7672
#	loc_2           -0.5652     -0.0376     -0.6241
#	loc_3            0.3793     -0.3232      0.1787
#	loc_4            0.7391     -0.5682      0.5909
#	loc_5           -0.5652     -0.0376     -0.6241
#	loc_6            0.6327     -0.1575      0.5748
#	           All:  0.1847      0.0310      0.2099
#	-----------------------------------------------

compute_Fstats <- function (data,infile=NA,MAF_threshold=0.0) {   # compute the multi-locus estimate of F_ST for diploid data, following Weir (1996)
  
  if (!is.na(infile)) data <- read.table(infile)                  # read the data
  rownames(data) <- NULL																				  # rename the rows
  colnames(data) <- c("ind","pop",paste("loc_",seq(1,ncol(data) - 2),sep = ""))	# rename the columns
  
  gen <- subset(data,select = -c(ind,pop))												# remove the first two columns (ind. ID and sample ID)
  pop <- subset(data,select = pop)																# take the vector of sample ID
  
  lst.pops <- unique(pop)																					# get the list of sample ID
  nbr.pops <- length(lst.pops)																		# compute the number of samples
  nbr.loci <- ncol(gen)																						# compute the number of loci
  counts <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci)						# define the matrix of allele counts
  nbr.hmzgtes.p.1 <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci) 	# define the matrix of type-1 homozygotes
  nbr.hmzgtes.p.2 <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci) 	# define the matrix of type-2 homozygotes
  n <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci)									# define the matrix of sample sizes
  
  cpt <- 1																										
  for (i in lst.pops) {																						          # loop over samples
    counts[cpt,] <- colSums(as.matrix(gen[which(pop == i),]),na.rm=TRUE)		# compute the allele counts from the dataset (per sample and per locus)
    nbr.hmzgtes.p.1[cpt,] <- colSums(as.matrix(gen[which(pop == i),] == 0),na.rm=TRUE)	# compute the number of type-1 homozygotes
    nbr.hmzgtes.p.2[cpt,] <- colSums(as.matrix(gen[which(pop == i),] == 2),na.rm=TRUE)	# compute the number of type-2 homozygotes
    n[cpt,] <- colSums(as.matrix(gen[which(pop == i),] != 9),na.rm=TRUE)								# compute the NUMBER OF INDIVIDUALS (per sample and per locus)
    cpt <- cpt + 1
  }
  
  n. <- colSums(n)																							  # compute the total NUMBER OF INDIVIDUALS (per locus)
  n2 <- colSums(n^2)																							# compute the sum of squared sample sizes (per locus)
  r <- nrow(unique(pop))																					# compute the number of sampled demes
  nc <- (n. - n2 / n.) / (r - 1.0)																# compute the n_c term in Weir (1996)
  
  p.1 <- (2 * n - counts) / (2 * n)																# compute the allele frequency for type-1 alleles (per sample and per locus)
  p.2 <- counts / (2 * n)																					# compute the allele frequency for type-2 alleles (per sample and per locus)
  
  
  vec.p.1.bar <- colSums(2 * n - counts) / colSums(2 * n)					# compute the overall allele frequency for type-1 alleles (per locus)
  vec.p.2.bar <- colSums(counts) / colSums(2 * n)									# compute the overall allele frequency for type-2 alleles (per locus)
  
  p.1.bar <- replicate(r,vec.p.1.bar)													# this is required to perform matrix operations afterwards
  p.2.bar <- replicate(r,vec.p.2.bar)													# this is required to perform matrix operations afterwards
  if (nbr.loci>1){
    p.1.bar <- t(p.1.bar)													# this is required to perform matrix operations afterwards
    p.2.bar <- t(p.2.bar)													# this is required to perform matrix operations afterwards
  }
  
  
  frq.hmzgtes.p.1 <- nbr.hmzgtes.p.1 / n													# compute the frequency of homozygotes of type 1
  frq.hmzgtes.p.2 <- nbr.hmzgtes.p.2 / n													# compute the frequency of homozygotes of type 1
  
  SSG <- colSums(n * (p.1 - frq.hmzgtes.p.1)) + colSums(n * (p.2 - frq.hmzgtes.p.2)) # compute the sum of squares for genes within individuals (Table 5.4, p. 177 in Weir 1996)
  SSI <- colSums(n * (p.1 + frq.hmzgtes.p.1 - 2 * p.1^2)) + colSums(n * (p.2 + frq.hmzgtes.p.2 - 2 * p.2^2)) # compute the sum of squares for individuals within populations (Table 5.4, p. 177 in Weir 1996)
  SSP <- 2 * colSums(n * (p.1 - p.1.bar)^2) + 2 * colSums(n * (p.2 - p.2.bar)^2) # compute the sum of squares for populations (Table 5.4, p. 177 in Weir 1996)
  
  MSG <- SSG / n.																								  # compute the observed mean squares for genes within individuals (Table 5.4, p. 177 in Weir 1996)
  MSI <- SSI / (n. - r)																						# compute the observed mean squares for individuals within populations (Table 5.4, p. 177 in Weir 1996)
  MSP <- SSP / (r - 1.0)																					# compute the observed mean squares for populations (Table 5.4, p. 177 in Weir 1996)
  
  F_ST_locus <-        (MSP - MSI) / (MSP + (nc - 1) * MSI + nc * MSG)   # compute the F_ST at each locus 
  F_IS_locus <-        (MSI - MSG) / (MSI + MSG)                         # compute the F_IS at each locus
  #F_IT_locus <- 1 - (2 * nc * MSG) / (MSP + (nc - 1) * MSI + nc * MSG)   # compute the F_IT at each locus
  

  maf <- maf_filter(p.1,MAF_threshold)
  
  
  F_ST_numerator   <- (MSP - MSI) 
  F_ST_denominator <- (MSP + (nc - 1) * MSI + nc * MSG)
    
  F_ST <- sum(MSP[maf$test] - MSI[maf$test]) / sum(MSP[maf$test] + (nc[maf$test] - 1) * MSI[maf$test] + nc[maf$test] * MSG[maf$test])		# compute the multilocus F_ST (see Weir 1996, p. 178)
  F_IS <- sum(MSI[maf$test] - MSG[maf$test]) / sum(MSI[maf$test] + MSG[maf$test])                         # compute the multilocus F_IS

  
  if (length(lst.pops)==2){
    F_C.numerator   <- (p.1[1,] - p.1[2,])^2
    F_C.denominator <- ((p.1[1,] + p.1[2,]) / 2 - p.1[1,] * p.1[2,])
    F_C_locus       <- F_C.numerator/F_C.denominator
    F_C             <- sum(F_C.numerator[maf$test])/sum(F_C.denominator[maf$test])
    if (nbr.loci==1){
      Fstats <- list(F_ST = F_ST_locus,
                     F_IS = F_IS_locus,
                     F_C  = F_C_locus,
                     numerator_F_ST   = (MSP - MSI),
                     denominator_F_ST = (MSP + (nc - 1) * MSI + nc * MSG),
                     numerator_F_IS   = (MSI - MSG),
                     denominator_F_IS = (MSI + MSG),
                     p = p.1,
                     n = n,
                     maf=maf)
    }else{
      Fstats <- list(F_ST = F_ST,
                     F_IS = F_IS,
                     F_C  = F_C,
                     F_ST_numerator = F_ST_numerator,
                     F_ST_denominator = F_ST_denominator,
                     F_ST_locus = F_ST_locus,
                     F_IS_locus = F_IS_locus,
                     F_C_locus  = F_C_locus,
                     p = p.1,
                     n = n,
                     maf=maf)
    }
  }else{
    Fstats <- list(F_ST = F_ST,
                   F_IS = F_IS,
                   F_ST_locus = F_ST_locus,
                   F_IS_locus = F_IS_locus,
                   F_ST_numerator = F_ST_numerator,
                   F_ST_denominator = F_ST_denominator,
                   p = p.1,
                   n = n,
                   maf=maf)
  }
  return (Fstats)
}


compute_freq <- function (data,infile=NA) {   # compute allele frequencies
  
  if (!is.na(infile)) data <- read.table(infile)                  # read the data
  rownames(data) <- NULL																				  # rename the rows
  colnames(data) <- c("ind","pop",paste("loc_",seq(1,ncol(data) - 2),sep = ""))	# rename the columns
  
  gen <- subset(data,select = -c(ind,pop))												# remove the first two columns (ind. ID and sample ID)
  pop <- subset(data,select = pop)																# take the vector of sample ID
  
  lst.pops <- unique(pop)																					# get the list of sample ID
  nbr.pops <- length(lst.pops)																		# compute the number of samples
  nbr.loci <- ncol(gen)																						# compute the number of loci
  counts <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci)						# define the matrix of allele counts
  n <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci)									# define the matrix of sample sizes
  
  cpt <- 1																										
  for (i in lst.pops) {																						          # loop over samples
    counts[cpt,] <- colSums(as.matrix(gen[which(pop == i),]),na.rm=TRUE)		# compute the allele counts from the dataset (per sample and per locus)
    n[cpt,] <- colSums(as.matrix(gen[which(pop == i),] != 9),na.rm=TRUE)		# compute the NUMBER OF INDIVIDUALS (per sample and per locus)
    cpt <- cpt + 1
  }
  
  p <- (2 * n - counts) / (2 * n)																# compute the allele frequency for type-1 alleles (per sample and per locus)

  return (p)
}















compute_FST_from_counts <- function(counts,multiallelic = TRUE) {
  r <- ncol(counts) / 2
  l <- seq(1,(2 * r),2)
  ss <- counts[,l] + counts[,(l + 1)]	
  ss2 <- rowSums((counts[,l] + counts[,(l + 1)])^2)
  n <- rowSums(ss)
  nc <- (n - ss2 / n) / (r - 1.0);
  p <- counts[,l] / ss
  q <- counts[,(l + 1)] / ss
  pbar <- rowSums(counts[,l]) / rowSums(ss)
  qbar <- rowSums(counts[,(l + 1)]) / rowSums(ss) 
  SSI <- rowSums(ss * (p - p^2) + ss * (q - q^2))
  SSP <- rowSums(ss * (p - pbar)^2 + ss * (q - qbar)^2)		
  MSI <- SSI / (n - r);
  MSP <- SSP / (r - 1.0);
  if (multiallelic) {
    F_ST <- sum(MSP - MSI)  / sum(MSP + (nc - 1) * MSI)
  }
  else {
    F_ST <- (MSP - MSI)  / (MSP + (nc - 1) * MSI)
  }
  return(F_ST)
}


compute_FST_from_genotype_counts <- function(genotypes) {
  
  
  # genotypes<-rbind(c(23,5,22),c(18,2,30))
  
  nbr.pops <- 2	# number of samples
  nbr.loci <- 1 # number of loci
  counts <- c(genotypes[1,1]*2+genotypes[1,2],
              genotypes[2,1]*2+genotypes[2,2])	# allele counts
  nbr.hmzgtes.p.1 <- genotypes[,1] 	# number of type-1 homozygotes
  nbr.hmzgtes.p.2 <- genotypes[,3] 	# number of type-2 homozygotes
  n <- rowSums(genotypes) 					# sample sizes (NUMBER OF INDIVIDUALS)
  
  n. <- sum(n) # compute the total NUMBER OF INDIVIDUALS (per locus)
  n2 <- sum(n^2) # compute the sum of squared sample sizes (per locus)
  r <- 2				# number of sampled demes
  nc <- (n. - n2 / n.) / (r - 1.0)		# compute the n_c term in Weir (1996)
  
  p.1 <- (2 * n - counts) / (2 * n)	# compute the allele frequency for type-1 alleles (per sample and per locus)
  p.2 <- counts / (2 * n)	 # compute the allele frequency for type-2 alleles (per sample and per locus)
  
  
  p.1.bar <- sum(2 * n - counts) / sum(2 * n)					# compute the overall allele frequency for type-1 alleles (per locus)
  p.2.bar <- sum(counts) / sum(2 * n)									# compute the overall allele frequency for type-2 alleles (per locus)
  
  frq.hmzgtes.p.1 <- nbr.hmzgtes.p.1 / n													# compute the frequency of homozygotes of type 1
  frq.hmzgtes.p.2 <- nbr.hmzgtes.p.2 / n													# compute the frequency of homozygotes of type 1
  
  SSG <- sum(n * (p.1 - frq.hmzgtes.p.1)) + sum(n * (p.2 - frq.hmzgtes.p.2)) # compute the sum of squares for genes within individuals (Table 5.4, p. 177 in Weir 1996)
  SSI <- sum(n * (p.1 + frq.hmzgtes.p.1 - 2 * p.1^2)) + sum(n * (p.2 + frq.hmzgtes.p.2 - 2 * p.2^2)) # compute the sum of squares for individuals within populations (Table 5.4, p. 177 in Weir 1996)
  SSP <- 2 * sum(n * (p.1 - p.1.bar)^2) + 2 * sum(n * (p.2 - p.2.bar)^2) # compute the sum of squares for populations (Table 5.4, p. 177 in Weir 1996)
  
  MSG <- SSG / n.																								  # compute the observed mean squares for genes within individuals (Table 5.4, p. 177 in Weir 1996)
  MSI <- SSI / (n. - r)																						# compute the observed mean squares for individuals within populations (Table 5.4, p. 177 in Weir 1996)
  MSP <- SSP / (r - 1.0)																					# compute the observed mean squares for populations (Table 5.4, p. 177 in Weir 1996)
  
  F_ST <-        (MSP - MSI) / (MSP + (nc - 1) * MSI + nc * MSG)   # compute the F_ST at each locus 
  F_IS <-        (MSI - MSG) / (MSI + MSG)                         # compute the F_IS at each locus

  return(F_ST)
}



prob.genotypes <- function(p,Fis){
  q   <- 1 - p
  PAA <- p^2 + Fis * p * q
  PAB <- 2*p*q*(1-Fis)
  PBB <- q^2 + Fis * p * q
  return(cbind(PAA,PAB,PBB))
}

