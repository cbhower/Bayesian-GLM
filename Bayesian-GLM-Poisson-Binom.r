
sink("hw5.txt")
cat("
    model {
    
    # Priors
    beta0 ~ dnorm(0, .00001)
    beta1 ~ dnorm(0, .00001)  

    # Likelihood    for (i in 1:n){
    counts[i] ~ dpois(mu[i])          
    log(mu[i]) <- beta0 +  beta1*PctCover[i]
    }
    # Deviance Observed
    for (i in 1:n){
    LL[i] <- -1*mu[i] + counts[i]*log(mu[i]) -
    logfact(counts[i])
    }
    # Deviance ideal
    for (i in 1:n){
    CountPred[i] ~ dpois(mu[i])
    LLP[i] <- -1*mu[i] + CountPred[i]*log(mu[i]) -
    logfact(CountPred[i])
    }
    # Monitoring
    dev <- -2*sum(LL[])
    devP <- -2*sum(LLP[])
    test <- step(dev-devP)
    bpvalue <- mean(test)
    #derived parameter
    bt.beta <- exp(.01*beta1)
    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(n = as.numeric(length(sal$Count)),
                 counts = as.numeric(sal$Count),
                 PctCover = as.numeric(sal$PctCover))

# Initial values
inits <- function() list(beta0 = runif(1, -2, 2), beta1 = runif(1, -2, 2))
# Parameters monitored
params <- c("beta0", "beta1", "dev", "devP", "bpvalue", "bt.beta")

# MCMC settings
ni <- 102000 
nt <- 50 
nb <- 20000 
nc <- 3

# Call WinBUGS from R 
out <- bugs(win.data, inits, params, "hw5.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
            bugs.directory = bugs.dir, working.directory = getwd())

1B)
# Autocorrelation, deviance plot, summary
acf(out$sims.matrix)
print(out , dig = 3)
dev.obs <- out$sims.matrix[,3]
dev.fit <- out$sims.matrix[,4]
plot(dev.fit ~ dev.obs , ylim=c(100, 240))
abline(a = 0 , b = 1 , col="red" , lwd = 3 )

1C)
sink("hw5.2.txt")
cat("
    model {
    # Priors
    beta0 ~ dnorm(0, .00001)
    beta1 ~ dnorm(0, .00001)  
    sd.alpha ~ dunif(0, 10)
    site.prec <- 1/(sd.alpha*sd.alpha)
    # Likelihood
    for (i in 1:n){
    counts[i] ~ dpois(mu[i])          
    log(mu[i]) <- beta0 +  beta1*PctCover[i] + alpha[S[i]]
    }
    for (i in 1:n){
    alpha[i] ~ dnorm(0 , site.prec)
    }
    # Deviance for dataset
    for (i in 1:n){
    LL[i] <- -1*mu[i] + counts[i]*log(mu[i]) -
    logfact(counts[i])
    }
    # Deviance for ideal datasets
    for (i in 1:n){
    CountPred[i] ~ dpois(mu[i])
    LLP[i] <- -1*mu[i] + CountPred[i]*log(mu[i]) -
    logfact(CountPred[i])
    }
    # Objects for Bayesian P value
    dev <- -2*sum(LL[])
    devP <- -2*sum(LLP[])
    test <- step(dev-devP)
    bpvalue <- mean(test)
    #derived parameter
    bt.beta <- exp(.01*beta1)
    }
    ",fill = TRUE)
sink(
# Bundle data
win.data <- list(n = as.numeric(length(sal$Count)),
                 counts = as.numeric(sal$Count),
                 PctCover = as.numeric(sal$PctCover),
                 S = as.numeric(sal$Site))
# Initial values
inits <- function() list(beta0 = runif(1, -2, 2), beta1 = runif(1, -2, 2), alpha = rnorm(length(sal$Site), 0 , 2))
# Parameters monitored
params <- c("beta0", "beta1", "dev", "devP", "bpvalue",
            "bt.beta" , "sd.alpha")

# MCMC settings
ni <- 102000 
nt <- 50 
nb <- 20000 
nc <- 3
# Call WinBUGS from R 
out2 <- bugs(win.data, inits, params, "hw5.2.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
             bugs.directory = bugs.dir, working.directory = getwd())


alf$N <- alf$Tumors + alf$Notumor
sink("hw5.3.txt")
cat("
    Model{
    #priors
    beta0~dnorm(0,0.00001)
    beta1~dnorm(0,0.0001)
    #likelihood
    for(i in 1:n){
    Tumors[i]~dbin(p[i],N[i])
    logit(p[i])<-beta0+beta1*lDose[i]
    }
    # Deviance Observed
    for(i in 1:n){
    LL[i]<- logfact(N[i]) - logfact(Tumors[i]) - 
    logfact(N[i] - Tumors[i]) + 
    Tumors[i]*log(p[i]) + 
    (N[i] - Tumors[i])*log(1-p[i])
    }
    # Deviance Fit
    for(i in 1:n){
    TumorsPred[i]~dbin(p[i],N[i])
    LLP[i]<- logfact(N[i]) - logfact(TumorsPred[i]) - 
    logfact(N[i] - TumorsPred[i]) +
    TumorsPred[i]*log(p[i]) +
    (N[i] - TumorsPred[i])*log(1-p[i])
    }
    dev<- -2*sum(LL[])
    devP<- -2*sum(LLP[])
    test<- step(dev - devP)
    bpvalue<-mean(test)
    }
    ", fill=TRUE)
sink()


# Bundle data
win.data <- list(Tumors = as.numeric(alf$Tumors),
                 n = as.numeric(length(alf$Tumors)),
                 N = as.numeric(alf$N),
                 lDose = as.numeric(alf$lDose))
# Initial values
inits <- function() list(beta0 = runif(1, -2, 2), beta1 = runif(1, -2, 2))
# Parameters monitored
params <- c("beta0", "beta1","dev","devP", "bpvalue")

# MCMC settings
ni <- 80000
nt <- 30
nb <- 20000
nc <- 3
# Call WinBUGS
out3 <- bugs(win.data, inits, params, "hw5.3.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, 
             bugs.directory = bugs.dir, working.directory = getwd())
2B)
# ACF, Deviance Plot, Summary
acf(out3$sims.matrix)
dev.obs3 <- out3$sims.matrix[,3]
dev.fit3 <- out3$sims.matrix[,4]
plot(dev.fit3 ~ dev.obs3 , xlim = c(75,180) , ylim = c(75,180))
abline(a = 0 , b = 1 , col="red" , lwd = 3 )
print(out3 , dig = 3)


2C)
sink("hw5.4.txt")
cat("
    Model{
    #priors
    beta0~dnorm(0,0.00001)
    beta1~dnorm(0,0.0001)
    sd.alpha ~ dunif(0 , 10)
    alpha.prec <- 1/(sd.alpha*sd.alpha)
    
    for (i in 1:n){
    alpha[i] ~ dnorm(0 , alpha.prec)
    }
    #likelihood
    for(i in 1:n){
    Tumors[i]~dbin(p[i],N[i])
    logit(p[i])<-beta0+beta1*lDose[i]+ alpha[T[i]]
    }
    # Observed Deviance
    for(i in 1:n){
    LL[i]<- logfact(N[i]) - logfact(Tumors[i]) - 
    logfact(N[i] - Tumors[i]) + 
    Tumors[i]*log(p[i]) + 
    (N[i] - Tumors[i])*log(1-p[i])
    }
    # Ideal Deviance
    for(i in 1:n){
    TumorsPred[i]~dbin(p[i],N[i])
    LLP[i]<- logfact(N[i]) - logfact(TumorsPred[i]) - 
    logfact(N[i] - TumorsPred[i]) +
    TumorsPred[i]*log(p[i]) +
    (N[i] - TumorsPred[i])*log(1-p[i])
    }
    dev<- -2*sum(LL[])
    devP<- -2*sum(LLP[])
    test<- step(dev - devP)
    bpvalue<-mean(test)
    bt.beta <- exp(beta1)
    }
    ", fill=TRUE)
sink()
# Bundle data
win.data <- list(Tumors = as.numeric(alf$Tumors),
                 n = as.numeric(length(alf$Tumors)),
                 N = as.numeric(alf$N),
                 lDose = as.numeric(alf$lDose),
                 T = as.numeric(alf$Tank))


# Initial values
inits <- function() list(beta0 = runif(1, -2, 2), beta1 = runif(1, -2, 2), alpha = rnorm(length(alf$Tank), 0 , 2) )
# Parameters monitored
params <- c("beta0", "beta1","dev","devP", "bpvalue", "sd.alpha","bt.beta")
# MCMC settings
ni <- 220000
nt <- 100
nb <- 20000
nc <- 3
# Call WinBUGS 
out4 <- bugs(win.data, inits, params, "hw5.4.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, 
             bugs.directory = bugs.dir, working.directory = getwd())

