library(reshape2)
library(fastDummies)
library(scales)

DYR <- read.csv("DYR_ECOLI.PF00186.scoremat.csv")
DYR_dis <- read.table("DYR_ECOLI-1ddr.dist.txt") 
DYR_dis <- DYR_dis[order(DYR_dis[,1]),]
DYR_energy <- read.table("DYR_ECOLI.fasta.mt.tsv") 
wild <- read.table("DYR_wild.txt",sep=",")

amino <- apply(data, 2, function(x){sort(unique(x))})
site <- sort(unique(DYR_energy[,1]))+1
site_number <- dim(data)[2]
l <- sapply(amino,length)
l <- c(0,l)
X <- dummy_cols(data[,1])[,-1]
for (i in 2:site_number) {
  X <- cbind(X,dummy_cols(data[,i])[,-1])
  print(i)
}
X <- as.matrix(X)
test <- unique(DYR_energy[,1:2])
test <- test[order(test[,1]),]
wild_type <- test[,2]
w <- diag(0,site_number)
for (i in 1:(site_number-1)) {
  for (j in (i+1):site_number) {
    a <- which(DYR_dis[,1]==(i -1)& DYR_dis[,2]==(j-1))
    w[i,j] <- DYR_dis[a,3] -> w[j,i]
  }
}

# create new data for CV
new_data_generate <- function(data,site){
  new_data <- data
  for (i in 1:site) {
    t <- table(new_data[,i])
    t <- t[which(t <5)]
    i_a <- as.numeric(names(t))
    for (j in 1:length(i_a)) {
      r_n <- rep(which(data[,i]==i_a[j]),ceiling(5/t[j])-1)
      new_data <- rbind(new_data,data[r_n,])
    }
  }
  return(as.matrix(new_data))
}
new_data <- new_data_generate(data,site_number)
new_X <- X_matrix(new_data,site_number)
new_X <- dummy_cols(new_data[,1])[,-1]
for (i in 2:site_number) {
  new_X <- cbind(new_X,dummy_cols(new_data[,i])[,-1])
}
X <- as.matrix(X)
new_X <- as.matrix(new_X)

