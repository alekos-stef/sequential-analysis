# Stein’s method
set.seed(1)
n <- 1000;n0 <- 10;mu <- 1;sigma <- 3;a <- 0.05;j <- 1
t <- qt(1-a/2,df=n0-1);mat.Stein <- matrix(,6,5)
Ts <- matrix() # matrix with the sample sizes T

for(d in c(1, 0.75 ,0.5, 0.25, 0.15, 0.05)) { # radius of c.i.
  s<-0;k <- qnorm(1-a/2,mean=0,sd=1)^2*sigma^2/d^2
  for(i in 1:n) {
      X1 <- rnorm(n0, mean=mu, sd=sigma)
      Ts[i] <- max(n0, floor(t^2*sd(X1)^2/d^2)+1)
      if(Ts[i]>n0) {
        X2 <- rnorm(Ts[i]-n0,mean=mu,sd=sigma) # choice of observ.
        X3 <- c(X1,X2) # matrix with the final sample
      }
      else {
        X3 <- X1
      }
      if(mu>=mean(X3)-d & mu<=mean(X3)+d) {
        s <- s+1
      }
    }
    mat.Stein [j,] <- c(k, mean(Ts), sd(Ts), mean(Ts)/k, s/n)
    j <- j+1
  }
print(mat.Stein)