# ASN for the model AR(1) for three values of ρ
set.seed(1);alpha <- 0.05;beta <- 0.1;theta0 <- 1;theta1 <- 1.4;sigma <-2
a<- log(beta/(1-alpha));b<- log((1-beta)/alpha);n <- 1000
mat.ASN.p <- matrix(,nrow=length(seq(1,1.4,0.04)),ncol=4,byrow=TRUE)
mat.ASN.p[,1] <-seq(1,1.4,0.04)
mat.OC.p <- matrix(,nrow=length(seq(1,1.4,0.04)),ncol=3,byrow=TRUE);k <-2

for(p in c(0.1,0.5,0.9)) {
  
  j <-1;ASN <- c();OC <- c()
  for(theta in seq(1,1.4,0.04)) {
    sigmaE <- sqrt((1-p^2))*sigma; oc <- c();count <- c()
    
    for(i in 1:n) {
      oc[i] <- 0
      En <- rnorm(1,0,sigmaE)
      count[i] <- 2
      S <- (theta1-theta0)/((1+p)*sigma^2)*(((1-p)*theta+En) - (1-p)*(theta1+theta0)/2 )
      
        while (S>a & S<b) {
          En <- rnorm(1,0,sigmaE)
          S <- S+(theta1-theta0)/((1+p)*sigma^2)*(((1-p)*theta+En) - (1-p)*(theta1+theta0)/2 )
          count[i] <- count[i]+1
        }
      
      if(S<a) {
        oc[i] <- oc[i]+1
      }
    }
    ASN[j] <- mean(count)
    OC[j] <- mean(oc)
    j <- j+1
  }
  mat.ASN.p[,k] <- ASN;mat.OC.p[,k-1] <- OC
  k <- k+1
}
totalmatrix<-cbind(mat.ASN.p,mat.OC.p);print(totalmatrix) # Final Table