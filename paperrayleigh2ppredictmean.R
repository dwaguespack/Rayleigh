ray2p.sample <- function(n,a,b){
  return(sqrt(-2*(b^2)*log(1-runif(n,0,1))) + a)
}

mle.a.b <- function(data){
  n <- length(data)
  xbar <- mean(data)
  h <- function(a){return((2*n^2)*sum(xbar-a)/sum((data-a)^2) - sum((data-a)^(-1)))}
  mlea <- uniroot(h,interval=c(min(data)-12*sqrt(2/(4-pi))*sd(data)/sqrt(n), min(data)),lower = min(data)-12*sqrt(2/(4-pi))*sd(data)/sqrt(n), upper = min(data))$root
  mleb <- sqrt(sum((data-mlea)^2)/(2*n))
  return(c(mlea,mleb))
}

std.mle <- function(){
  mat.a <- mat.b <- matrix(rep(0,10000*50),10000,50)
  for(i in 5:50){
    for(j in 1:10000){
      sam <- ray2p.sample(i,0,1)
      mles <- mle.a.b(sam)
      mat.a[j,i] <- mles[1]
      mat.b[j,i] <- mles[2]
    }
    print(i)
  }
  write.csv(mat.a,file="ray2p.mle.a.emalg.10k",row.names = F)
  write.csv(mat.b,file="ray2p.mle.b.emalg.10k",row.names = F)
}

zbar.dist <- function(){
  dist.zbar <- matrix(0,10000,20)
  for(j in 1:20){
    for(i in 1:10000){
      dist.zbar[i,j] <- mean(ray2p.sample(j,0,1))
    }
  }
  return(dist.zbar)
}


pred.int.fmean.mle <- function(data,m,cl){
  n <- length(data)
  xbar <- mean(data)
  s <- sd(data)
  mles <- mle.a.b(data)
  ah <- mles[1]
  bh <- mles[2]
  Qa <- ah - bh*ahs[,n]/bhs[,n]
  Qb <- bh/bhs[,n]
  Qxb.m <- Qa + Qb*zbar[,m]
  return(quantile(Qxb.m,c(.5*(1-cl),.5*(1+cl))))
}

cov.fmean.mle <- function(a,b,n,m,cl,sims){
  cov <- 0
  
  for(i in 1:sims){
    fmean <- mean(ray2p.sample(m,a,b))
    pi <- pred.int.fmean.mle(ray2p.sample(n,a,b),m,cl)
    if(pi[1]<fmean & fmean<pi[2]){cov <- cov+1}
  }
  return(cov/sims)
}

cov.plot <- function(a,n,m,cl,sims){
  b = seq(from=0.1, to=3, by=0.1)
  cov <- rep(0,length(b))
  for(i in 1:length(b)){
    cov[i] <- cov.fmean.mle(a,b[i],n,m,cl,sims)
    print(i)
  }
  plot(b,cov,ylim = c(.5,1),type = "l")
  lines(b,rep(.95,length(b)),type="l",lty=2)
}







std.ray.mme <- function(){
  A <- matrix(rep(0,10000),10000,50)
  B <- matrix(rep(0,10000),10000,50)
  for(i in 5:50){
    for(j in 1:10000){
      sam <- ray2p.sample(i,0,1)
      xbar <- mean(sam)
      s <- sd(sam)
      A[j,i] <- xbar - s*sqrt(pi/(4-pi))
      B[j,i] <- sqrt(2/(4-pi))*s
    }
    print(i)
  }
  write.csv(A,"std.ray.mme.a.10k",row.names = F)
  write.csv(B,"std.ray.mme.b.10k",row.names = F)
}



pred.int.fmean.mme <- function(data,m,cl){
  n <- length(data)
  xbar <- mean(data)
  s <- sd(data) 
  mles <- mle.a.b(data)
  ah <- xbar - s*sqrt(pi/(4-pi))
  bh <- sqrt(2/(4-pi))*s
  Qa <- ah - bh*ahs[,n]/bhs[,n]
  Qb <- bh/bhs[,n]
  Qxb.m <- Qa + Qb*zbar[,m]
  return(quantile(Qxb.m,c(.5*(1-cl),.5*(1+cl))))
}

cov.fmean.mme <- function(a,b,n,m,cl,sims){
  cov <- 0
  
  for(i in 1:sims){
    fmean <- mean(ray2p.sample(m,a,b))
    pi <- pred.int.fmean.mme(ray2p.sample(n,a,b),m,cl)
    if(pi[1]<fmean & fmean<pi[2]){cov <- cov+1}
  }
  return(cov/sims)
}



L.est.std <- function(){
  est.mat.a <- matrix(0,10000,50)
  est.mat.b <- matrix(0,10000,50)
  for(i in 5:50){
    for(j in 1:10000){
      sam <- ray2p.sample(i,0,1)
      sam <- sort(sam)
      l1 <- mean(sam)
      l2 <- 2*sum(((1:i)-1)*sam)/(i*(i-1)) - l1
      sqrt2 <- sqrt(2)
      ah <- l1-sqrt2*l2/(sqrt2-1)
      lamh <- (gamma(3/2)^2)*(3-2*sqrt2)/(2*l2^2)
      bh <- 1/sqrt(2*lamh)
      est.mat.a[j,i] <- ah
      est.mat.b[j,i] <- bh
    }
    print(i)
  }
  write.csv(est.mat.a, file="std.ray.lmom.a",row.names=F)
  write.csv(est.mat.b, file="std.ray.lmom.b",row.names=F)
}




pred.int.fmean.lme <- function(data,m,cl){
  n <- length(data)
  data <- sort(data)
  l1 <- mean(data)
  l2 <- 2*sum(((1:n)-1)*data)/(n*(n-1)) - l1
  sqrt2 <- sqrt(2)
  ah <- l1-sqrt2*l2/(sqrt2-1)
  lamh <- (gamma(3/2)^2)*(3-2*sqrt2)/(2*l2^2)
  bh <- 1/sqrt(2*lamh)
  Qa <- ah - bh*ahs[,n]/bhs[,n]
  Qb <- bh/bhs[,n]
  Qxb.m <- Qa + Qb*zbar[,m]
  return(quantile(Qxb.m,c(.5*(1-cl),.5*(1+cl))))
}


cov.fmean.lme <- function(a,b,n,m,cl,sims){
  cov <- 0
  
  for(i in 1:sims){
    fmean <- mean(ray2p.sample(m,a,b))
    pi <- pred.int.fmean.lme(ray2p.sample(n,a,b),m,cl)
    if(pi[1]<fmean & fmean<pi[2]){cov <- cov+1}
  }
  return(cov/sims)
}









