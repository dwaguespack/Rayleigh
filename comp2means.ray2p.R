ray2p.sample <- function(n,a,b){
  return(sqrt(-2*(b^2)*log(1-runif(n,0,1))) + a)
}

mle.a.b <- function(data){
  n <- length(data)
  xbar <- mean(data)
  h <- function(a){return((2*n^2)*(xbar-a)/sum((data-a)^2) - sum((data-a)^(-1)))}
  mlea <- uniroot(h,interval=c(min(data)-12*sqrt(2/(4-pi))*sd(data)/sqrt(n), min(data)),lower = min(data)-12*sqrt(2/(4-pi))*sd(data)/sqrt(n), upper = min(data))$root
  mleb <- sqrt(sum((data-mlea)^2)/(2*n))
  return(c(mlea,mleb))
}

std.mle <- function(){
  mat.a <- mat.b <- matrix(rep(0,10000*30),10000,30)
  for(i in 5:30){
    for(j in 1:10000){
      sam <- ray2p.sample(i,0,1)
      mles <- mle.a.b(sam)
      mat.a[j,i] <- mles[1]
      mat.b[j,i] <- mles[2]
    }
    print(i)
  }
  write.csv(mat.a,file="ray2p.mle.a.10k",row.names = F)
  write.csv(mat.b,file="ray2p.mle.b.10k",row.names = F)
}

ci.ray2p.diffmu.mle <- function(data1,data2,alpha){
  n1 <- length(data1)
  n2 <- length(data2)
  mles1 <- mle.a.b(data1)
  mle.a1 <- mles1[1]
  mle.b1 <- mles1[2]
  mles2 <- mle.a.b(data2)
  mle.a2 <- mles2[1]
  mle.b2 <- mles2[2]
  Qa1 <- mle.a1 - mle.b1*ahs1[,n1]/bhs1[,n1]
  Qb1 <- mle.b1/bhs1[,n1]
  Qa2 <- mle.a2 - mle.b2*ahs2[,n2]/bhs2[,n2]
  Qb2 <- mle.b2/bhs2[,n2]
  Qmu1 <- Qa1 + Qb1*sqrt(pi/2)
  Qmu2 <- Qa2 + Qb2*sqrt(pi/2)
  Qdiffmu <- Qmu1-Qmu2
  return(quantile(Qdiffmu,c(alpha/2,1-alpha/2)))
}

cov.diffmu.mle <- function(a1,b1,a2,b2,n1,n2,alpha,sims){
  cov <- 0
  len <- rep(0,sims)
  diffmu <- a1+b1*sqrt(pi/2) - (a2+b2*sqrt(pi/2))
  for(i in 1:sims){
    data1 <- ray2p.sample(n1,a1,b1)
    data2 <- ray2p.sample(n2,a2,b2)
    ci <- ci.ray2p.diffmu.mle(data1,data2,.05)
    if(ci[1]<diffmu & diffmu<ci[2]){cov <- cov+1}
    len[i] <- ci[2]-ci[1]
  }
  return(c(round(cov/sims,3),round(mean(len),2)))
} 


####################################################################################################

std.ray.mme <- function(){
  A <- matrix(rep(0,10000),10000,30)
  B <- matrix(rep(0,10000),10000,30)
  for(i in 5:30){
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


ci.diff.mean.mme <- function(dataX,dataY,cl){
  n <- length(dataX); m <- length(dataY)
  xbar <- mean(dataX); ybar <- mean(dataY)
  sx <- sd(dataX); sy <- sd(dataY)
  phix <- sqrt(pi)+sqrt(2)*ahs1[,n]/bhs1[,n]
  phiy <- sqrt(pi)+sqrt(2)*ahs2[,m]/bhs2[,m]
  ax <- xbar-sx*phix/sqrt(4-pi)
  ay <- ybar-sy*phiy/sqrt(4-pi)
  bx <- sqrt(2/(4-pi))*sx/bhs1[,n]
  by <- sqrt(2/(4-pi))*sy/bhs2[,m]
  diff.mean <- ax-ay+sqrt(pi/2)*(bx-by)
  return(quantile(diff.mean,c(.5*(1-cl),.5*(1+cl))))
}

cov.ci.diff.mean <- function(aX,bX,aY,bY,n,m,cl,sim){
  mux <- aX+sqrt(pi/2)*bX
  muy <- aY+sqrt(pi/2)*bY
  diff <- mux-muy
  cov <- 0
  len <- rep(0,sim)
  for(i in 1:sim){
    ci <- ci.diff.mean.mme(ray2p.sample(n,aX,bX),ray2p.sample(m,aY,bY),cl)
    if(ci[1]<diff & diff<ci[2]){cov <- cov+1}
    len[i] <- ci[2] - ci[1]
  }
  return(c(round(cov/sim,3), round(mean(len),2)))
}