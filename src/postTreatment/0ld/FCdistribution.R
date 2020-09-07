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
Fc_selected <- list()

for (simulation in 1:nrow(sim_table) ){
  sim <- sim_table$simID[simulation]
  message(paste("Calculating footprint of selection for scenario",sim,"\n"))
  
  Fc_m2    <- array(NA,number_of_replicates)
  Fc_m1    <- numeric()
  Fc_m1_MAF<- numeric()
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      Fc_m2[replic] <- SNP_list$FC_obs[which(SNP_list$type=="m2")]
      Fc_m1         <- c(Fc_m1,SNP_list$FC_obs[which(SNP_list$type!="m2")])
      loci_to_sample<- intersect(which(SNP_list$type!="m2"),which(!is.na(SNP_list$FC_p_value)))
      Fc_m1_MAF     <- c(Fc_m1_MAF,SNP_list$FC_obs[loci_to_sample])
      
    }
  }
  
  Fc_neutral[[length(Fc_neutral)+1]] <- Fc_m1
  Fc_neutralMAF[[length(Fc_neutralMAF)+1]] <- Fc_m1_MAF
  Fc_selected[[length(Fc_selected)+1]] <- Fc_m2
  
  
}

save(sim_table,Fc_neutral,Fc_neutralMAF,Fc_selected,file="Results/Fc_distribution.RData")




#load(file="Results/Fc_distribution.RData")
source("percentile.boxplot.R")





## Plot FC in function of selfing

selection_mode <- "NM"
sel_coef                  <- 0.5
selection_period_duration <- 25
simulations2plot          <- intersect(which(sim_table$sel_coef==sel_coef),
                                       which(sim_table$selection_period_duration==selection_period_duration))
simulations2plot          <- intersect(simulations2plot,
                                       which(sim_table$selection_mode==selection_mode))
sigma_values              <- sim_table$sigma[simulations2plot]

pdf_file_name <- paste0("Results/FC_s",sel_coef,"_dT",selection_period_duration,"_",selection_mode,".pdf")
pdf(file=pdf_file_name)
#jpg_file_name <- paste0("Results/FC_s",sel_coef,"_dT",selection_period_duration,"_",selection_mode,".jpg")
#jpeg(file=jpg_file_name)


boxplot(  Fc_neutral[[ simulations2plot[1] ]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.03,
          names=sigma_values[1],
          ylab = "Fc",
          xlab = expression("selfing rate, "*sigma),
          ylim = c(0,2),
          xlim = c(min(sigma_values),max(sigma_values)),
          main = substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )
boxplot(  Fc_neutralMAF[[ simulations2plot[1] ]],
          at=sigma_values[1],
          outline=FALSE,
          axes = FALSE,
          boxwex=0.01,
          col="green",
          names=sigma_values[1],
          add=T)
points( array(sigma_values[1],length(Fc_selected[[ simulations2plot[1] ]])),
        Fc_selected[[ simulations2plot[1] ]],
        col="red",cex=0.3)

for (i in 2:length(sigma_values)){
  boxplot(  Fc_neutral[[ simulations2plot[i] ]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.03,
            names=sigma_values[i],
            add=T)
  boxplot(  Fc_neutralMAF[[ simulations2plot[i] ]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.01,
            col="green",
            names=sigma_values[i],
            add=T)
  points( array(sigma_values[i],length(Fc_selected[[ simulations2plot[i] ]])),
          Fc_selected[[ simulations2plot[i] ]],
          col="red",cex=0.3)
}
axis(side=1,at=sigma_values)
if(selection_mode=="NM") texto<-"New mutation"
if(selection_mode=="SV") texto<-"Standing variation"
text(0.5,1.5,adj=0.5,texto,cex=3,col="grey")

dev.off()
