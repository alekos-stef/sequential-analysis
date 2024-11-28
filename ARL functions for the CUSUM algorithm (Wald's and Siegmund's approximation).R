set.seed(1)
theta0 <-1;theta1<-2;sigma<-1;n<-1000
h<-4;l<-1;N<-seq(0.9,2.5,0.05)
ARLsim<-matrix();ARLwald<-matrix()
ARLsieg<-matrix();ARLmat<-matrix(,nrow=length(N),ncol=4);meanEZ<-matrix()

for(theta in N) {
  
  k<-matrix()
  for(j in 1:n) {
    
    S<-matrix();Z<-matrix()
    X<-rnorm(1,mean=theta,sd=sigma)
    S[1]<-0;Z[1]<-(theta1-theta0)/sigma^2*(X-(theta1+theta0)/2)
    S[2]<-S[1]+Z[1]
    i<-2
    
    while(S[i]<h) {
      
      if(S[i]<=0) {
        S[i]<-0;
      }
      
      X<-rnorm(1,mean=theta,sd=sigma)
      Z[i]<-(theta1-theta0)/sigma^2*(X-(theta1+theta0)/2)
      S[i+1]<-S[i]+Z[i]
      i<-i+1
    }
    k[j]<-i
  }
  
  ARLsim[l] <- mean(k)
  ARLwald[l]<-ifelse(theta==1.5,h^2,1/(theta-1.5)*(h+exp(-(2*theta-3)*h)/(2*theta-3)-1/(2*theta-3)))
  ARLsieg[l]<-ifelse(theta==1.5,(h+1.166)^2,1/(theta-1.5)*(h+1.166+exp(-(2*theta-3)*(h+1.166))/(2*theta-3)-1/(2*theta-3)))
  meanEZ[l] <-(theta1-theta0)/sigma^2*(theta-(theta1+theta0)/2)
  l<-l+1
}

ARLmat[,1]<-meanEZ;ARLmat[,2]<-ARLsim;ARLmat[,3]<-ARLwald;ARLmat[,4]<-ARLsieg;
par(mar=c(3,3,2,2))
plot(N,ARLsieg,type="l",lty=4,xlab="",xlim=c(0.9,2.1),ylab="",main="")
points(N,ARLsim,pch=4)
lines(N,ARLwald,type="l",lty=1)