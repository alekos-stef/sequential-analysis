# Purely Sequential Method
set.seed(1)
n <- 1000; n0 <- 10; mu <- 1; sigma <- 3; a <- 0.05; l <-1
Ta <- matrix() # matrix with T values
mat.Seq <- matrix(,6,5); # final table of algorithm
X <- matrix() # table with the final sample
for(d in c(1, 0.75 , 0.5, 0.25, 0.15, 0.05)) { 
  # radius of c.i.
  k <- qnorm(1-a/2,mean=0,sd=1)^2*sigma^2/d^2
  s <- 0
    for(j in 1:n) {
      i <- n0
      X <- rnorm(n0, mean=mu, sd=sigma)
      k0 <- qnorm(1-a/2, mean=0, sd=1)^2*sd(X)^2/d^2
      while(i<k0) { # loop for estimation of Ta
        i <- i+1
        X[i] <- rnorm(1, mean=mu, sd=sigma)
        k0 <- qnorm(1-a/2, mean=0, sd=1)^2*sd(X)^2/d^2
      }
      Ta[j] <- i
        if(mu>=mean(X)-d & mu<=mean(X)+d) {
          s <- s+1
        }
      }
    mat.Seq[l,] <- c(k, mean(Ta), sd(Ta), mean(Ta)/k, s/n)
    l <- l+1
  }
print(mat.Seq)