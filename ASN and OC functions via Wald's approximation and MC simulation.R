set.seed(1) # Wald’s approximation
alpha <- 0.05;beta <- 0.1;a <- log(beta/(1-alpha));b <- log((1-beta)/alpha)
theta0 <- 1;theta1 <- 1.4;sigma <- 2;k <- (theta1-theta0)/sigma^2;theta.hat <-(theta1+theta0)/2
i <- 1;A.wald <- matrix(nrow=length(seq(1,1.4,0.04)), ncol=4, byrow=TRUE)
Q <- c();ASN <-c();EZ1 <-c();t0<-c();EZ1.sq <- c()

for(theta in seq(1,1.4,0.04)) {
  
  if (theta!=theta.hat) {
    t0[i] <- 2/(theta1-theta0)*(theta-theta.hat)
    Q[i] <- (exp(1)^(-t0[i]*b)-1 )/(exp(1)^(-t0[i]*b) - exp(1)^(-t0[i]*a))
    EZ1[i] <- (theta1-theta0)/sigma^2*(theta-theta.hat)
    ASN[i] <- (a*Q[i]+b*(1-Q[i]))/EZ1[i]
  } else {
  t0[i]<-0
  Q[i] <- b/(b-a)
  EZ1.sq[i] <- k^2*(sigma^2+theta^2)+k^2*theta.hat*(theta.hat-2*theta)
  ASN[i] <- (a^2*Q[i]+b^2*(1-Q[i]))/EZ1.sq[i]
  }
  A.wald[i,] <- c(t0[i],theta,Q[i], ASN[i])
  i <- i+1
}
j<-1
ASN.sim <- c();oc.sim <-c(); # ASN and OC simulation

  for(theta.sim in seq(1,1.4,0.04)) {
    
    count <- c();q <- c();n <- 1000
    
    for(i in 1:n) {
      
      count[i] <- 1;q[i] <- 0
      X <- rnorm(1,mean=theta.sim, sd=sigma)
      Sx <- log(dnorm(X,mean=theta1,sd=sigma)/dnorm(X, mean=theta0, sd=sigma))
        while(Sx>a & Sx<b) {
          X1 <- rnorm(1, mean=theta.sim, sd=sigma)
          Sx <- Sx +log(dnorm(X1, mean=theta1, sd=sigma)/dnorm(X1, mean=theta0, sd=sigma))
          count[i] <- count[i]+1
        }
      if(Sx<=a) {
        q[i] <- q[i]+1
      }
    }
  
    ASN.sim[j] <- mean(count)
    oc.sim[j] <- mean(q)
    j <- j+1
  }
A.sim<-cbind(oc.sim,ASN.sim); A<-cbind(A.wald,A.sim); print(A) # Final Table