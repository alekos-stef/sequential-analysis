# Srivastava’s Method
set.seed(1)
library(MASS)
p=2;m=10;Mean=c(1,2);S=matrix(c(1,0.5,0.5,2),nrow=2,ncol=2);
alpha=0.05;n=1000
u<-p*(m-1)/(m-p)*qf(1-alpha, df1=p, df2=m-p);mat.Sriva <- matrix(,ncol=5,nrow=6);j<-1
a<- qchisq(1-alpha, df=p);
lambda_S <- max(eigen(S)$values)

    for(d in c(1,0.75, 0.5, 0.25, 0.15, 0.05)) { # radius of region
      C <- qchisq(1-alpha, df=p)*lambda_S/d^2 # optimal sample
      sum=0;XTbar <- c();T <- c();X<-c()
      for(i in 1:n) {
        X <- mvrnorm(m, mu=Mean, Sigma=S)
        Sn <- cov(X) # sample covariance matrix
        lambda_n <- max(eigen(Sn)$values)
        N <- m # Αρχική τιμή του N
        while(N<a/d^2*lambda_n){
          N <- N+1
          Xnew <- mvrnorm(1, mu=Mean, Sigma=S)
          X <- rbind(X,Xnew)
          Sn <- cov(X)
          lambda_n <- max(eigen(Sn)$values)
        }
        T[i] <- N # estimation of sample size
        XTbar <- c(mean(X[,1]),mean(X[,2])) # mean vector
        if (t(XTbar-Mean)%*%(XTbar-Mean)<=d^2){ # confidence region
          sum=sum+1
        }
      }
      mat.Sriva[j,] <- c(C, mean(T), sd(T), mean(T)/C, sum/n)
      j <- j+1
    }
print(mat.Sriva)