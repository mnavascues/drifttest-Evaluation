#######################################################
#
# Evaluation of effective population size estimates (effect of MAF)
#
#######################################################


# read simulations description table
sim_table <- read.table("Simulations/Fc.simparams",header=T)

source("slim_tools.R")


number_of_replicates <- 100
N <- 500
sample_size <- 50
  
Proportion_fixed_sample1 <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
Proportion_fixed_sample2 <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
Proportion_fixed_sample1MAF <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
Proportion_fixed_sample2MAF <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))

for (simulation in 1:nrow(sim_table) ){
  sim <- sim_table$simID[simulation]
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      Proportion_fixed_sample1[replic,simulation] <- length(union(which(SNP_list$derived.x==100),which(SNP_list$ancestral.x==100)))/nrow(SNP_list)
      Proportion_fixed_sample2[replic,simulation] <- length(union(which(SNP_list$derived.y==100),which(SNP_list$ancestral.y==100)))/nrow(SNP_list)
      Proportion_fixed_sample1MAF[replic,simulation] <- length(union(which(SNP_list$derived.x[maf]==100),which(SNP_list$ancestral.x[maf]==100)))/length(which(maf))
      Proportion_fixed_sample2MAF[replic,simulation] <- length(union(which(SNP_list$derived.y[maf]==100),which(SNP_list$ancestral.y[maf]==100)))/length(which(maf))
    }
  }

}


save(sim_table,
     Proportion_fixed_sample1,
     Proportion_fixed_sample2,
     Proportion_fixed_sample1MAF,
     Proportion_fixed_sample2MAF,
     file="ProportionFixed.RData")















#load(file="ProportionFixed.RData")
source("percentile.boxplot.R")

# Make figure in function of selfing rate

sel_coef                  <- 0.005
selection_period_duration <- 25
simulations2plot <- intersect(which(sim_table$sel_coef==sel_coef),which(sim_table$selection_period_duration==selection_period_duration))
sigma_values              <- sim_table$sigma[simulations2plot]

pdf_file_name <- paste0("ProportionLociFixedSample1_s",sel_coef,"_dT",selection_period_duration,".pdf")
pdf(file=pdf_file_name)

y_limits <- c( 0,1 )

boxplot(  Proportion_fixed_sample1[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.02,
          names=sigma_values[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selfing rate, "*sigma),
          ylim = y_limits,
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
for (i in 2:length(sigma_values)){
  boxplot(  Proportion_fixed_sample1[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.02,
            names=sigma_values[i],
            add=T)

}
axis(side=1,at=sigma_values)

dev.off()



pdf_file_name <- paste0("ProportionLociFixedSample2_s",sel_coef,"_dT",selection_period_duration,".pdf")
pdf(file=pdf_file_name)

y_limits <- c( 0,1 )

boxplot(  Proportion_fixed_sample2[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.02,
          names=sigma_values[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selfing rate, "*sigma),
          ylim = y_limits,
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
for (i in 2:length(sigma_values)){
  boxplot(  Proportion_fixed_sample2[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.02,
            names=sigma_values[i],
            add=T)
  
}
axis(side=1,at=sigma_values)

dev.off()





pdf_file_name <- paste0("ProportionLociFixedSample1MAF_s",sel_coef,"_dT",selection_period_duration,".pdf")
pdf(file=pdf_file_name)

y_limits <- c( 0,1 )

boxplot(  Proportion_fixed_sample1MAF[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.02,
          names=sigma_values[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selfing rate, "*sigma),
          ylim = y_limits,
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
for (i in 2:length(sigma_values)){
  boxplot(  Proportion_fixed_sample1MAF[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.02,
            names=sigma_values[i],
            add=T)
  
}
axis(side=1,at=sigma_values)

dev.off()



pdf_file_name <- paste0("ProportionLociFixedSample2MAF_s",sel_coef,"_dT",selection_period_duration,".pdf")
pdf(file=pdf_file_name)

y_limits <- c( 0,1 )

boxplot(  Proportion_fixed_sample2MAF[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.02,
          names=sigma_values[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selfing rate, "*sigma),
          ylim = y_limits,
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
for (i in 2:length(sigma_values)){
  boxplot(  Proportion_fixed_sample2MAF[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.02,
            names=sigma_values[i],
            add=T)
  
}
axis(side=1,at=sigma_values)

dev.off()


