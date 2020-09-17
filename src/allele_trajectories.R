    source("src/fun/ColorBlindPalette.R")
    palette(cbPalette1)
    transp_color <- makeTransparent("black",alpha=0.2)
      
    traj <- function(p,sigma,w11,w12,w22,tgen=25){
      p.array<-p
      f <- sigma/(2*(1-sigma/2))  
      for(i in 1:tgen){
        x11 <- p^2 + f*(1-p)*p
        x12 <- 2*(1-p)*p*(1-f)
        x22 <- (1-p)^2 + f*(1-p)*p
        wbar<- w11*x11 +w12*x12 + w22 * x22
        p <- (w11*x11 + w12*x12)/wbar 
        p.array<-c(p.array,p)
      }
      return(p.array)
    }
      
      
    svg(file="results/AlleleTrajectory.svg",width=4.5,height=7)
    par(mar=c(4,4,0.2,1)+0.1)
    layout(matrix(1:2, 2,1,byrow = TRUE))
    
    load(paste0("results/scenario007/replicates/scenario007_1_data.RData"))
    plot(0:25,c(1/1000,trajectory/1000),
      ylim=c(0,1),
      xlab="generations",ylab="allele frequency",
      type="l",col=transp_color)
    for(i in 2:100){
      #print(i)
      load(paste0("results/scenario007/replicates/scenario007_",i,"_data.RData"))
      if (length(trajectory)<25){
        trajectory <- c(trajectory,rep(1000,25-length(trajectory)))
      }
      lines(0:25,c(1/1000,trajectory/1000),col=transp_color)
    }
    text(0,0.95,"a",cex=2)
      
    load(paste0("results/scenario031/replicates/scenario031_1_data.RData"))
    plot(0:25,c(0.1,1-trajectory/1000),
         ylim=c(0,1),
         xlab="generations",ylab="allele frequency",
         type="l",col=transp_color)
    for(i in 2:100){
      #print(i)
      load(paste0("results/scenario031/replicates/scenario031_",i,"_data.RData"))
      if (length(trajectory)<25){
        trajectory <- c(trajectory,rep(0,25-length(trajectory)))
      }
      lines(0:25,c(0.1,1-trajectory/1000),col=transp_color)
    }
    text(0,0.95,"b",cex=2)
      
    dev.off()
      
