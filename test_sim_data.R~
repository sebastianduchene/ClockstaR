library(ape)

tr.fix <- rtree(10)

rates1 <- abs(rnorm(18, sd = 0.1))
rates2 <- abs(rnorm(18, sd = 0.1))
rates3 <- abs(rnorm(18, sd = 0.1))
rates4 <- abs(rnorm(18, sd = 0.1))

trees.list <- list()

for(i in 1:20){
      trees.list[[i]] <- tr.fix
      if(i <= 5){
      	   trees.list[[i]]$edge.length <- abs(rates1 + rnorm(18, 0, 0.01)) 
      }else if(i > 5 && i <= 10){
      	    trees.list[[i]]$edge.length <- abs(rates2 + rnorm(18, 0, 0.01)) 
      }else if(i >= 10 && i < 15){
      	    trees.list[[i]]$edge.length <- abs(rates3 + rnorm(18, 0, 0.01))  
      }else{
	    trees.list[[i]]$edge.length <- abs(rates4 + rnorm(18, 0, 0.01)) 
      }
}

names(trees.list) <- paste0("tree", 1:20)
class(trees.list) <- "multiPhylo"

source("funs2/bsd.dist.R")
source("funs2/bsd.matrix.R")

library("doParallel")
library("foreach")
library(cluster)
source("funs2/bsd.matrix.para.R")

trees.bsd <- bsd.matrix(trees.list)

source("funs2/partitions.bsd.R")
