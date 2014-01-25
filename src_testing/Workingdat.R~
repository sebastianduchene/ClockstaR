getwd()
dir("./src/")

setwd("src")
for(i in dir()){
    if(!(any(grepl("~", i)))){
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

# Gap statistic fails to detect multiple clocks if the maximum number is 2

trees.2 <- trees.test[1:2]
di.2 <- min.dist.topo.mat(trees.2)
grs.2 <- get.all.groups(as.matrix(di.2[[1]]))
# Returns an error

trees.3 <- trees.test[c(1, 4, 9)]
di.3 <- min.dist.topo.mat(trees.3)
grs.3 <- get.all.groups(as.matrix(di.3[[1]]))
# This does not produce an erro, but it should.

trees.4 <- trees.test[c(1, 4, 9, 8)]
di.4 <- min.dist.topo.mat(trees.4)
grs.4 <- get.all.groups(as.matrix(di.4[[1]]))
# works!!
