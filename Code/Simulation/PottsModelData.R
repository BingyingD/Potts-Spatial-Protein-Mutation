library(fastDummies)
Gibbs.func <- function(z,Theta,Gamma,P,K){
  gap <- seq(0,by=K,length.out=P-1)
  for (i in 1:P) {
    Gamma_i <- matrix(Gamma[i,],nrow = K)
    X <- z[-i]+gap
    prob <- exp(rowSums(Gamma_i[,X])+Theta[i,]) 
    z[i] <- sample(K,1,prob=prob)
  }
  return(z)
}

W_generate.func <- function(P,sp){
  test <- matrix(rbeta(P^2,4,2),P)
  test <- (test+t(test))/2
  diag(test) <- 0
  test2 <- matrix(runif(sp^2,0,0.1),sp)
  test2 <- (test2+t(test2))/2
  diag(test2) <- 0
  test[1:sp,1:sp] <- test2
  return(test)
}

data_generate.func <- function(n,K,P,Gamma,R){
  # Gamma: P * K(P-1)K
  # initial
  z0 <- sample(K,P,T)
  result <- matrix(NA,R+1,P)
  result[1,] <- z0
  # updata z
  for (i in 1:R) {
    z0 <- Gibbs.func(result[i,],Theta,Gamma,P,K)
    result[i+1,] <- z0
  }
  result_burnin <- result[-(1:round(0.4*R)),]
  id <- sample(1:dim(result_burnin)[1],n)
  return(result_burnin[id,])
}

X_matrix <- function(data,site_number){
  X <- dummy_cols(data[,1])[,-1]
  for (i in 2:site_number) {
    X <- cbind(X,dummy_cols(data[,i])[,-1])
  }
  return(as.matrix(X))
}

n <- 1000
K <- 20
P <- 25
R <- 10000
sp <- 5
sk <- 5
Theta <- matrix(runif(P*K,0,1.5),P,K)
Theta[which(Theta<0.5)] <- 0
gamma_generate.func <- function(P,K,sp,sk){
  gamma <- cbind(replicate(sk,as.vector(cbind(replicate(sp-1,c(runif(sk,0.5,2),rep(0,K-sk))),
                                              matrix(0,K,P-sp)))),
                 matrix(0,(P-1)*K,K-sk))
  return(as.vector(t(gamma)))
}
Gamma <- t(cbind(replicate(sp,gamma_generate.func(P,K,sp,sk)),
                 matrix(0,(P-1)*K^2,P-sp)))
location <- apply(Gamma, 1, function(x){
  which(x!=0)
})
data <- data_generate.func(n,K,P,Gamma,R)
site_number <- P
amino <- apply(data, 2, function(x){sort(unique(x))})
l <- c(0,sapply(amino,length))
amino <- unlist(amino)
X <- X_matrix(data,site_number)
W <- W_generate.func(P,sp)
wild_type <- rep(1,P)

