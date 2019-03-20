ray2p.sample <- function(n,a,b){
  return(sqrt(-2*(b^2)*log(1-runif(n,0,1))) + a)
}

mles <- function(data,r){
  n <- length(data)
  ordX <- sort(data)
  cen.ordX <- ordX[(r+1):n]
  b <- function(a){return(sqrt(sum((cen.ordX-a)^2)/(2*(n-r)) + r*((ordX[r+1]-a)^2)/(2*(n-r))))}
  bsp <- function(a){return(-sum(cen.ordX-a)/(n-r) - r*(ordX[r+1]-a)/(n-r))}
  g <- function(a){return(sum((-b(a)^2)/(cen.ordX-a) + cen.ordX -a) + r*(ordX[r+1]-a))}
  gp <- function(a){return(sum((-bsp(a)*(cen.ordX-a)-b(a)^2)/(cen.ordX-a)^2 - 1) - r)}
  
  a0 <- mean(cen.ordX) - sd(cen.ordX)*sqrt(pi/(4-pi))
  a1 <- a0 - g(a0)/gp(a0)
  i <- 1
  repeat{
    if(abs(a1-a0)<1e-5 | i>50){break}
    a0 <- a1
    a1 <- a0-g(a0)/gp(a0)
    i <- i+1
  }
  ah <- a1
  bh <- b(ah)
  return(c(ah,bh))
}

std.ray.mles <- function(){
  mat.a <- matrix(rep(0,1000000*4),1000000,4)
  mat.b <- matrix(rep(0,1000000*4),1000000,4)
  r <- c(3,5,7,9)
  for(j in 1:4){
    for(i in 1:1000000){
      m <- mles(ray2p.sample(10,0,1),0)
      mat.a[i,j] <- m[1]
      mat.b[i,j] <- m[2]
    }
    print(j)
  }
  return(cbind(mat.a,mat.b))
}

ci.mean.cen <- function(data,r){
  
  mle <- mles(data,r)
  ah <- mle[1]
  bh <- mle[2]
  ahs <- std.mle.n10.r1[,1]
  bhs <- std.mle.n10.r1[,2]
  t <- 1
  Qa <- ah-bh*std.mle.n20.r0[,1]/std.mle.n20.r0[,2]
  Qb <- bh/std.mle.n20.r0[,2]
  Qmu <- Qa + Qb*sqrt(pi/2)
  #print(c(ah,bh))
  cv <- (quantile(Qmu,c(.05,.95)) - ah)/bh
  #print(quantile((sqrt(pi/2)-ahs)/bhs,c(.05,.95)))
  # print(quantile(Qmu,c(alpha/2,1-alpha/2)))
  # print(c(ah,bh))
  print(cv)
  #return(quantile(Qmu,c(.05,.95)))
}



cov.prob <- function(a,b){
  mu <- a+b*sqrt(pi/2)
  cov <- 0
  for(i in 1:10000){
    ci <- ci.mean.cen(ray2p.sample(10,a,b),7)
    if(ci[1]<mu & mu < ci[2]){cov <- cov+1}
  }
  return(cov/10000)
}

ray2p.nr <- function(data,iter){
  xbar <- mean(data)
  n <- length(data)
  a0 <- xbar - sqrt(((n-1)/n)*var(data))*sqrt(pi/(4-pi))
  gp <- function(a){
    return(((xbar-a)*2*n^2)/sum((data-a)^2) - sum((data-a)^(-1)))
  }
  gpp <- function(a){
    return(((-2*n^2)*sum((data-a)^2) + 2*n*(xbar-a))/(sum((data-a)^2)^2) - sum((data-a)^(-2)))
  }
  a1 <- a0 - gp(a0)/gpp(a0)
  i <- 1
  while(abs(a1-a0)>1e-7 & i < iter){
    a0 <- a1
    a1 <- a0 - gp(a0)/gpp(a0)
    i <- i+1
  }
  bh <- sqrt(1/(2*n/sum((data-a1)^2)))
  return(c(a1,bh))
}

mle.mu <- function(data,cl){
  alpha <- 1-cl
  n <- length(data)
  s <- sqrt((n-1)*var(data)/n)
  mles <- ray2p.nr(data,50)
  ah <- mles[1]
  bh <- mles[2]
  Qa <- ah - bh*ahsl[,n]/bhsl[,n]
  Qb <- bh/bhsl[,n]
  Qmu <- Qb*sqrt(pi/2)+Qa
  return((quantile(Qmu,c(alpha/2,1-alpha/2))-ah)/bh)
}
