#######################################################
#
# Calculating False Discovery Rate
#
#######################################################


# read simulations description table
sim_table <- read.table("Simulations/Fc.simparams",header=T)

#pdf(file="Figures.pdf",paper="a4r")

number_of_replicates <- 100
threshold_power <- 0.05

FDR <- array(NA,nrow(sim_table))

for (simulation in 1:nrow(sim_table) ){
  sim <- sim_table$simID[simulation]
  message(paste("Calculating footprint of selection for scenario",sim,"\n"))
  
  q_values_m1    <- numeric()
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      q_values_m1 <- c(q_values_m1,SNP_list$FC_q_value[intersect(which(SNP_list$type!="m2"),which(is.na(SNP_list$distance_bp)))])
    }
  }
  q_values_m1 <- q_values_m1[which(!is.na(q_values_m1))]
  FDR[simulation] <- length(which(q_values_m1<threshold_power))/length(q_values_m1)
  
}

sim_table <- cbind(sim_table,FDR)

save(sim_table,file="Results/FalseDiscoveryRate.RData")
#load(file="Results/FalseDiscoveryRate.RData")






#color blind palette
require(graphics)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette(cbPalette)

selection_mode <- "NM"
selection_period_duration <- 25
simulations2plot_total <- intersect(which(sim_table$selection_period_duration==selection_period_duration),
                                    which(sim_table$selection_mode==selection_mode))
sel_coef               <- unique(sim_table$sel_coef[simulations2plot_total])

#pdf_file_name <- paste0("Results/FDR_dT",selection_period_duration,"_",selection_mode,".pdf")
#pdf(file=pdf_file_name)
jpg_file_name <- paste0("Results/FDR_dT",selection_period_duration,"_",selection_mode,".jpg")
jpeg(file=jpg_file_name)



simulations2plot <- intersect(simulations2plot_total,which(sim_table$sel_coef==sel_coef[1]))
sigma_values     <- sim_table$sigma[simulations2plot]
FDR            <- sim_table$FDR[simulations2plot]
plot(sigma_values,FDR,
     type="l",
     col=1,
     ylim=c(0,0.1),
     ylab = expression("False Discovery Rate, "*italic(FDR)),
     xlab = expression("Selfing rate, "*sigma),
     main=substitute( italic(N)*"="*500*", "*italic(t)*"="*dt  ,list(dt=selection_period_duration)) ,
     lwd=3)
for (i in 2:length(sel_coef)){
  simulations2plot <- intersect(simulations2plot_total,which(sim_table$sel_coef==sel_coef[i]))
  sigma_values     <- sim_table$sigma[simulations2plot]
  FDR            <- sim_table$FDR[simulations2plot]
  lines(sigma_values,FDR,col=i,lwd=3)
}
if(selection_mode=="NM") texto<-"New mutation"
if(selection_mode=="SV") texto<-"Standing variation"
text(0.5,0.09,adj=0.5,texto,cex=3,col="grey")
dev.off()


