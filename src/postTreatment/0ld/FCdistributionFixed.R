#######################################################
#
# Fc distribution (neutral vs. selected)
#
#######################################################


# read simulations description table
sim_table <- read.table("Simulations/Fc.simparams",header=T)

#pdf(file="Figures.pdf",paper="a4r")

number_of_replicates <- 100

# read simulations description table
sim_table <- read.table("Simulations/Fc.simparams",header=T)


Fc_neutral  <- list()
Fc_neutralMAF  <- list()
Fc_neutralFixed  <- list()
Fc_selected <- list()

for (simulation in 1:nrow(sim_table) ){
  sim <- sim_table$simID[simulation]

  Fc_m2      <- array(NA,number_of_replicates)
  Fc_m1      <- numeric()
  Fc_m1_MAF  <- numeric()
  Fc_m1_fixed<- numeric()
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      Fc_m2[replic]   <- SNP_list$FC_obs[which(SNP_list$type=="m2")]
      Fc_m1           <- c(Fc_m1,SNP_list$FC_obs[which(SNP_list$type!="m2")])
      loci_to_sample  <- intersect(which(SNP_list$type!="m2"),which(maf))
      Fc_m1_MAF       <- c(Fc_m1_MAF,SNP_list$FC_obs[loci_to_sample])
      loci_to_sample  <- intersect(loci_to_sample,which(SNP_list$ancestral.x==100))
      Fc_m1_fixed     <- c(Fc_m1_fixed,SNP_list$FC_obs[loci_to_sample])
      
    }
  }
  
  Fc_neutral[[length(Fc_neutral)+1]] <- Fc_m1
  Fc_neutralMAF[[length(Fc_neutralMAF)+1]] <- Fc_m1_MAF
  Fc_selected[[length(Fc_selected)+1]] <- Fc_m2
  Fc_neutralFixed[[length(Fc_neutralFixed)+1]] <- Fc_m1_fixed
  
  
}

save(sim_table,Fc_neutral,Fc_neutralMAF,Fc_selected,Fc_neutralFixed,file="Fc_distributionFixed.RData")










## Plot FC in function of selfing

sel_coef                  <- 0.005
selection_period_duration <- 25
simulations2plot          <- intersect(which(sim_table$sel_coef==sel_coef),which(sim_table$selection_period_duration==selection_period_duration))
sigma_values              <- sim_table$sigma[simulations2plot]

pdf_file_name <- paste0("FCfixed_s",sel_coef,"_dT",selection_period_duration,".pdf")
pdf(file=pdf_file_name)


boxplot(  Fc_neutral[[ simulations2plot[1] ]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.03,
          names=sigma_values[1],
          ylab = "Fc",
          xlab = expression("selfing rate, "*sigma),
          ylim = c(0,0.5),
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
points( array(sigma_values[1],length(unique(Fc_neutralFixed[[ simulations2plot[1] ]]))),
        unique(Fc_neutralFixed[[ simulations2plot[1] ]]),
        col="red",cex=0.3)

for (i in 2:length(sigma_values)){
  boxplot(  Fc_neutral[[ simulations2plot[i] ]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.03,
            names=sigma_values[i],
            add=T)
  points( array(sigma_values[i],length(unique(Fc_neutralFixed[[ simulations2plot[i] ]]))),
          unique(Fc_neutralFixed[[ simulations2plot[i] ]]),
          col="red",cex=0.3)
}
axis(side=1,at=sigma_values)



dev.off()
