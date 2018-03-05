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
  
Ne_estimates <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
Ne_estimatesNNM <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
Fc_multilocus <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
Fc_multilocusNNM <- matrix(NA,nrow=number_of_replicates,ncol=nrow(sim_table))
expected_Ne <- array(NA,nrow(sim_table))

for (simulation in 1:nrow(sim_table) ){
  sim <- sim_table$simID[simulation]
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      Ne_estimates[replic,simulation]     <- Fc$Ne_FC_multi
      Fc_multilocus[replic,simulation]    <- Fc$FC_multi
      loci2sample <- intersect(union(which(SNP_list$derived.x==0),which(SNP_list$ancestral.x==0)),which(maf))
      Fc_multilocusNNM[replic,simulation] <- sum(SNP_list$FC_num[loci2sample])/sum(SNP_list$FC_denom[loci2sample])
      Ne_estimatesNNM[replic,simulation]  <- Compute.F_c.N_e(mean_FC=Fc_multilocusNNM[replic,simulation],
                                                             dT=sim_table$selection_period_duration[simulation],
                                                             S1=sample_size,
                                                             S2=sample_size)
      
      
    }
  }

  expected_Ne[simulation] <- (2-sim_table$sigma[simulation])*N
}

sim_table <- cbind(sim_table,expected_Ne)

save(sim_table,Ne_estimates, Fc_multilocus,Ne_estimatesNNM,Fc_multilocusNNM,file="EffectivePopSizeNNM.RData")















#load(file="EffectivePopSizeNNM.RData")
source("percentile.boxplot.R")


# Make figure in function of selfing rate

sel_coef                  <- 0.005
selection_period_duration <- 100
sigma_value               <- 0
simulation2plot          <-  intersect(which(sim_table$sigma==sigma_value),
                                        intersect(which(sim_table$sel_coef==sel_coef), 
                                                  which(sim_table$selection_period_duration==selection_period_duration)))

pdf_file_name <- paste0("EffectivePopSizeNNM_s",sel_coef,"_dT",selection_period_duration,"_sigma",sigma_value,".pdf")
pdf(file=pdf_file_name)

y_limits <- c( 800,1200 )

boxplot(  Ne_estimates[,simulation2plot],Ne_estimatesNNM[,simulation2plot],
          outline=FALSE,
          names= c("all loci","NNM filter"),
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = "",
          ylim = y_limits,
          main = paste0("N=500, s=",sel_coef,", t=",selection_period_duration,", sigma=",sigma_value) )
          #main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt", "*sigma*"="*sig ) ,list(ss=sel_coef,dt=selection_period_duration,sig=sigma_value)) )

abline(h=sim_table$expected_Ne[simulation2plot],col="red")
dev.off()




















# Make figure in function of selfing rate

sel_coef                  <- 0.005
selection_period_duration <- 100
simulations2plot <- intersect(which(sim_table$sel_coef==sel_coef),which(sim_table$selection_period_duration==selection_period_duration))
sigma_values              <- sim_table$sigma[simulations2plot]

pdf_file_name <- paste0("EffectivePopSizeNNM_s",sel_coef,"_dT",selection_period_duration,".pdf")
pdf(file=pdf_file_name)

y_limits <- c( 500,1500 )

boxplot(  Ne_estimates[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.02,
          names=sigma_values[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selfing rate, "*sigma),
          ylim = y_limits,
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
boxplot(  Ne_estimatesNNM[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          axes = FALSE,
          boxwex=0.01,
          add=T,
          col="green")

for (i in 2:length(sigma_values)){
  boxplot(  Ne_estimates[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.02,
            names=sigma_values[i],
            add=T)
  boxplot(  Ne_estimatesNNM[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.01,
            add=T,
            col="green")
}
lines(sigma_values,sim_table$expected_Ne[simulations2plot],col="red")
axis(side=1,at=sigma_values)
dev.off()




