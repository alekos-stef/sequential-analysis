# Healy’s Method
set.seed(1);library(MASS)
p=2;m=10;Mean=c(1,2);S=matrix(c(1,0.5,0.5,2),nrow=2,ncol=2);alpha=0.05;n=1000
u<-p*(m-1)/(m-p)*qf(1-alpha, df1=p, df2=m-p);
mat.Healy <- matrix(,ncol=5,nrow=6);j<-1
lambda_S <- max(eigen(S)$values) # eigenvalue of Σ matrix

    for(d in c(1,0.75, 0.5, 0.25, 0.15, 0.05)) { # radius of region
      C <- qchisq(1-alpha, df=p)*lambda_S/d^2
      sum=0;XTbar <- c();T <- c();X3<-c()
      for(i in 1:n) {
  
          X1 <- mvrnorm(m, mu=Mean, Sigma=S)
          Sn <- cov(X1) # sample covariance matrix
          lambda <- max(eigen(Sn)$values)
          T[i] <- max(m,floor(u*lambda/d^2)+1)
         
            if(T[i]>m) {
              X2 <- mvrnorm(T[i]-m, mu=Mean, Sigma=S)
              X3 <- rbind(X1,X2)
            } else {
                X3 <- X1
              }
          XTbar <- c(mean(X3[,1]),mean(X3[,2])) # mean value vector
          
          if (t(XTbar-Mean)%*%(XTbar-Mean)<=d^2){ # Confidence region
            sum=sum+1
          }
        }
      mat.Healy[j,] <- c(C, mean(T), sd(T), mean(T)/C, sum/n)
      j <- j+1
    }
print(mat.Healy)