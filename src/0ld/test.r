library(coda)
require(gtools)

true_Fis    <- 0.99
true_p      <- 0.2
sample_size <- 50

X <- rmultinom(1, sample_size, c( true_p*true_p+true_Fis*true_p*(1-true_p),
                                  2*(1-true_p)*true_p*(1-true_Fis),
                                  (1-true_p)*(1-true_p)+true_Fis*true_p*(1-true_p)) )
(X)
X_AA <- X[1,1]
X_Aa <- X[2,1]
X_aa <- X[3,1]

observed_p <- (2*X_AA+X_Aa)/(2*sample_size)

Fis  <- true_Fis

posterior <- function(p){
  gamma_AA <- p*p+Fis*p*(1-p)
  gamma_Aa <- 2*(1-p)*p*(1-Fis)
  gamma_aa <- (1-p)*(1-p)+Fis*p*(1-p)
  likelihood <- log(gamma_AA^X_AA) + log(gamma_Aa^X_Aa) + log(gamma_aa^X_aa)
  prior      <- dbeta(p,1,1,log=T)
  posterior  <- likelihood + prior
  posterior
}


#posterior(observed_p)

proposal <- function(p){
  p <- rbeta(1,p*sample_size/10,(1-p)*sample_size/10)
  p
}

proposal <- function(p){
  p <- rbeta(1,p*sample_size*2+1,(1-p)*sample_size*2+1)
  p
}


#proposal(observed_p)

MH_MCMC <- function(startvalue, iterations){
  chain    <- array(NA,iterations+1)
  chain[1] <- startvalue
  for (i in 1:iterations){
    prop  <- proposal(chain[i])
    proba <- exp( posterior(prop) - posterior(chain[i]) )
    if (runif(1) < proba){
      chain[i+1] <- prop
    }else{
      chain[i+1] <- chain[i]
    }
  }
  return(chain)
}

chain <- MH_MCMC(observed_p,10000)
plot(chain,type="l")

#acceptance rate
1-mean(duplicated(chain))

length(chain)
ESS <- effectiveSize(as.mcmc(chain))

(ESS)
ESS/length(chain)

chain_subsample <- chain[seq(1,length(chain),10)]
plot(chain_subsample,type="l")

den_chain <- density(chain)
#den_subchain <- density(chain_subsample)

plot(den_chain$x,den_chain$y,
     xlab="allele frequency", ylab = "probability density",
     xlim=c(0,1),ylim=c(0,12),type="l",lwd=2)
x <- seq(0,1,0.005)
den <- dbeta(x,shape1=2*X_AA+X_Aa+1,shape2=2*X_aa+X_Aa+1)
lines(x,den,type = "l",col="red",lty=2,lwd=2)
#lines(den_subchain$x,den_subchain$y,type = "l",col="blue",lty=3,lwd=2)



gamma <- rdirichlet(1000,c(X_AA+1,X_Aa+1,X_aa+1))
den_diricchlet <- density(apply(gamma, 1, function(x){x[1]+x[2]/2} ))
lines(den_diricchlet$x,den_diricchlet$y,type = "l",col="green",lty=3,lwd=2)

