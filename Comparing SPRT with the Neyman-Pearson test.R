set.seed(1)
sigma<-2;theta0<-1;theta1<-1.4;alpha<-0.05;n<-200;beta<-matrix();j<-1;N<-10000;

for(n in c(30,50,100,150,200,250,300)) {
  beta[j]<-pnorm((theta0-theta1)/sqrt(sigma^2/n)+qnorm(1-alpha,mean=0,sd=1),mean=0,sd=1)
  j<-j+1
}

matrixNP<-matrix(ncol=3,nrow=length(beta)) # final table N-P
matrixNP[,1]<-c(30,50,100,150,200,250,300)
matrixNP[,2]<-rep(alpha,length(beta));matrixNP[,3]<-beta
matrixSPRT<-matrix(nrow=length(beta),ncol=8)

for(j in 1:length(beta)) {
  A <- beta[j]/(1-alpha); B <- (1-beta[j])/alpha
  a<-log(A);b<-log(B)
  mat.H0<-matrix();mat.H1<-matrix() # tables with T0 and T1
  k<-0;l<-0 # counters for Pr(I) and Pr(II)
  
    for(i in 1:N) {
      X0<-rnorm(1,mean=theta0,sd=sigma)
      Z0<-(theta1-theta0)/sigma^2*(X0-(theta0+theta1)/2)
      Sn.0<-Z0 # initial value of sum Sn under H0
      countH0<-1 # counter for mean value E[T0]
      X1<-rnorm(1,mean=theta1,sd=sigma)
      Z1<-(theta1-theta0)/sigma^2*(X1-(theta0+theta1)/2)
      Sn.1<-Z1 # initial value of Sn under H1
      countH1<-1 # counter for the mean value E[T1]
      
      while((a<Sn.0) & (Sn.0<b)) { # loop for estimation Pr[I] and E[T0]
        X0<-rnorm(1,mean=theta0,sd=sigma)
        Z0<-(theta1-theta0)/sigma^2*(X0-(theta0+theta1)/2)
        Sn.0<-Sn.0+Z0
        countH0<-countH0+1
      }
      
      while((a<Sn.1) & (Sn.1<b)) { # loop for estim. of Pr[II] and E[T1]
        X1<-rnorm(1,mean=theta1,sd=sigma)
        Z1<-(theta1-theta0)/sigma^2*(X1-(theta0+theta1)/2)
        Sn.1<-Sn.1+Z1
        countH1<-countH1+1
      }
      
      if(Sn.0>=b) {
        k<-k+1
      }
      if(Sn.1<=a) {
        l<-l+1
      }
      
     mat.H0[i]<-countH0
     mat.H1[i]<-countH1
     
    }
  
    matrixSPRT[j,]<-c(mean(mat.H0),sd(mat.H0)/length(mat.H0)
                      ,mean(mat.H1),sd(mat.H1)/length(mat.H1),k/N,l/N,k/N+l/N,alpha+beta[j])
   }
 totalmatrix<-cbind(matrixNP,matrixSPRT) # final table
 colnames(totalmatrix)<-c("n","alpha","beta_NP","Mean[T;H0]", "s.e(T0_bar)","Mean[T;H1]", "s.e(T1_bar)","Pr[I]","Pr[II]", "Pr[I]+Pr[II]","a+b")
 print(totalmatrix)
 par(mfrow=c(1,2));par(mar=c(3,3,2,2))
hist(mat.H0,breaks=40,xlab="",ylab="",main="");hist(mat.H1,breaks=40,xlab="",
                                                      ylab="",main="")