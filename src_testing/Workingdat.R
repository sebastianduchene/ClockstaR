getwd()
dir("./src/")

setwd("src")
for(i in dir()){
    if(any(grepl("~", i))){
        source(i)
    }
 }
setwd("..")



###
# get.all.groups should print an error when the data are not in the correct format

library(phangorn)
library(cluster)

trees.test <- read.tree("testing_trees.trees")

tree.dists <- min.dist.topo.mat(trees.test)

grs1 <- get.all.groups(tree.dists[[1]], pam.results = T)


grs2 <- get.all.groups(as.matrix(tree.dists[[1]]), pam.results = T)

source("src/get.all.groups.R")

######### SOLVED
