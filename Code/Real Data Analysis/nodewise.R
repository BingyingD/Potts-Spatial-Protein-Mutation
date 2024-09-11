# install.packages("remotes") # uncomment if remotes is not already installed
# remotes::install_github("nielsrhansen/sglOptim", build_vignettes = TRUE)
# remotes::install_github("nielsrhansen/msgl", build_vignettes = TRUE)

library(msgl) #load the package

args <- commandArgs(trailingOnly = TRUE)
ii1 <- as.numeric(args[[1]]) # nodewise regression

CV_grid_func <- function(new_x,new_classes,group,gWeights,
                         alpha_seq,lambda1,fold,sWeights_new){
  result <- sapply(alpha_seq, function(x){
    cv1 <- msgl::cv(new_x, new_classes, grouping = group,groupWeights =gWeights,
                    alpha = x,standardize = F,sampleWeights = sWeights_new,
                    lambda=lambda1,fold=fold)
    Error <- Err(cv1)
    newl <-lambda1[order(Error)[1]]
    error <- Error[order(Error)[1]]
    return(c(newl,error))
  })
  ind <- which.min(result[2,])
  return(data.frame(newl=result[1,ind],alpha=alpha_seq[ind]))
}

load("DYR.RData") 
sWeights_new <- read.csv("sWeights_new.txt") # calculate sample weights from .m file 
sWeights <- read.csv("sWeights.txt") 

i <- site[ii1]
x <- X[,-((cumsum(l)[i]+1):cumsum(l)[i+1])]
group <- as.factor(rep(1:(site_number-1),times=l[-c(1,i+1)]))
classes <- as.factor(data[,i])
new_classes <-as.factor(new_data[,i]) #data & the relative classes for CV

W <- w[-i,i] # distance between sites 
sig <- sd(W)
number_elements <- as.numeric(table(group))
n <- dim(data)[1]
p <- length(number_elements)

# Our method
# group weights with distances
gWeights2 <-  sqrt(number_elements/n)+sqrt(2*log(p)/n)
gWeights <- ( 1-exp(-W^2/(sig^2))) * gWeights2 
new_x <-  new_X[,-((cumsum(l)[i]+1):cumsum(l)[i+1])]
# grid search for CV
lambda1 <- seq(2,0.001,length.out=20)
cv1 <- CV_grid_func(new_x,new_classes,group,gWeights,
                                seq(0,1,.1),lambda1,fold=5,sWeights_new)
result <- fit(x,classes,grouping = group,groupWeights =gWeights,
              standardize = F,lambda =c(cv1$newl*100,cv1$newl),
              alpha = cv1$alpha,sampleWeights = sWeights)$beta
beta_hat <- matrix(sapply(result, as.matrix)[,2], nrow = l[site[ii1]+1])

write.csv(beta_hat,file=paste0("/Archive/DYR/OurMethod_V",i,".csv"),row.names = F)

# SGL
cv2 <- CV_grid_func(new_x,new_classes,group,gWeights2,
                    seq(0,1,.1),lambda1,fold=5,sWeights_new)
result2 <- fit(x,classes,grouping = group,groupWeights =gWeights2,
              standardize = F,lambda =c(cv2$newl*100,cv2$newl),
              alpha = cv2$alpha,sampleWeights = sWeights)$beta
beta_hat_SGL <- matrix(sapply(result2, as.matrix)[,2], nrow = l[site[ii1]+1])
write.csv(beta_hat_SGL,
          file=paste0("/Archive/DYR/SGL_V",i,".csv"),row.names = F)





