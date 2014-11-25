\name{bsd.matrix.para}
\alias{bsd.matrix.para}
\title{
Estimate the sBSDmin distance for all pairs of trees in a list with parallelisation. 
}
\description{
bsd.matrix estimates the sBSDmin distance for all paris of trees in a list. It returns an object of class bsd with the pairswise distances and the scaling factors. The scaling factors represent the differences in the magnitude of the rates among trees. See the help file for bsd.dist for more details.
}
\usage{
bsd.matrix.para(tree.list, para = F, ncore = 1)
}
\arguments{
  \item{tree.list}{
A list of trees. It can be a list where each item is a tree, or an object of class multiPhylo. The trees can have names. See the example bellow.
}
  \item{para}{
logical. select T to run in parallel.
}
  \item{ncore}{
number of cores for parallelisation.
}
}

\value{
An object of class 'bsd', which is a list with 
\item{comp1 }{An object of class dist with the pairwise distances for trees}
\item{comp2 }{An object of class matrix with the scaling factors among trees}
}

\author{
Sebastian Duchene
}

\seealso{
bsd.dist
bsd.matrix
}
\examples{

# Create a list of trees of class multiPhylo with four patterns of among-lineage rate variation 
tr.fix <- rtree(10)

set.seed(12345)
rates1 <- abs(rnorm(18, sd = 0.1))
set.seed(123456)
rates2 <- abs(rnorm(18, sd = 0.1))
set.seed(1234567)
rates3 <- abs(rnorm(18, sd = 0.1))
set.seed(12345678)
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

\dontrun{
# Estimate sBSDmin distance for all pairs of trees:
trees.bsd <- bsd.matrix.para(trees.list, para = T, ncore = 2)
}


## The function is currently defined as
function (tree.list, para = F, ncore = 1) 
{
require(foreach)
require(doParallel)

    if (length(tree.list) <= 3) {
        stop("The number of gene trees is <= 3. ClockstaR requires at least gene 4 trees")
    }
    bsd.dist <- function(tree1, tree2) {
        list.tr <- list()
        list.tr[[1]] <- tree1
        list.tr[[2]] <- tree2
        lens <- c(sum(tree1$edge.length), sum(tree2$edge.length))
        tree1 <- list.tr[lens == max(lens)][[1]]
        tree2 <- list.tr[lens == min(lens)][[1]]
        tree.dist.opt <- function(x) {
            tree3 <- tree2
            tree3$edge.length <- tree2$edge.length * x
            return(dist.topo(tree1, tree3, method = "score"))
        }
        opt.dist <- optim(0, fn = tree.dist.opt, method = "Brent", 
            lower = 0, upper = 50)
        min.bdi <- opt.dist$value
        scaling <- opt.dist$par
        tree2.scaled <- tree2
        tree2.scaled$edge.length <- tree2$edge.length * scaling
        root.scaling <- 0.05/mean(c(tree1$edge.length[tree1$edge.length > 
            1e-05], tree2.scaled$edge.length[tree2.scaled$edge.length > 
            1e-05]))
        tree1.root.scaled <- tree1
        tree2.root.scaled <- tree2.scaled
        tree1.root.scaled$edge.length <- tree1$edge.length * 
            root.scaling
        tree2.root.scaled$edge.length <- tree2.scaled$edge.length * 
            root.scaling
        min.bdi.root.scaled <- dist.topo(tree1.root.scaled, tree2.root.scaled, 
            method = "score")
        res.vect <- c(min.bdi.root.scaled, scaling, min.bdi)
        names(res.vect) <- c("min.bdi.scaled", "scaling.factor", 
            "min.bdi")
        return(res.vect)
    }
    sub.trees <- list()
    for (k in 2:length(tree.list)) {
        sub.trees[[k]] <- tree.list[1:k - 1]
    }
    compute.tree.dists <- function(tree.sub.list, fix.tree) {
        res <- sapply(tree.sub.list, function(a) {
            return(bsd.dist(fix.tree, a))
        })
        return(res)
    }
    if (para == T) {
        cl <- makeCluster(ncore)
        registerDoParallel(cl)
        print(paste("Clusters registered as follows: ", cl))
        res.par <- foreach(s.trees = sub.trees, j = 1:length(tree.list), 
            .packages = "ape") \%dopar\% compute.tree.dists(tree.sub.list = s.trees, 
            fix.tree = tree.list[[j]])
        stopCluster(cl)
    }
    else if (para == F) {
        res.par <- foreach(s.trees = sub.trees, j = 1:length(tree.list), 
            .packages = "ape") \%do\% compute.tree.dists(tree.sub.list = s.trees, 
            tree.list[[j]])
    }
    res.list <- list()
    res.list[[1]] <- matrix(NA, nrow = length(tree.list), ncol = length(tree.list))
    for (m in 2:nrow(res.list[[1]])) {
        res.list[[1]][m, 1:ncol(res.par[[m]])] <- res.par[[m]][1, 
            ]
    }
    rownames(res.list[[1]]) <- names(tree.list)
    colnames(res.list[[1]]) <- names(tree.list)
    res.list[[1]] <- as.dist(res.list[[1]])
    res.list[[2]] <- matrix(NA, nrow = length(tree.list), ncol = length(tree.list))
    for (m in 2:nrow(res.list[[2]])) {
        res.list[[2]][m, 1:ncol(res.par[[m]])] <- res.par[[m]][2, 
            ]
    }
    rownames(res.list[[2]]) <- names(tree.list)
    colnames(res.list[[2]]) <- names(tree.list)
    class(res.list) <- "bsd"
    return(res.list)
  }
}
\keyword{ sBSDmin }
