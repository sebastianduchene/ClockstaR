\name{min.dist.topo.mat.para}
\alias{min.dist.topo.mat.para}
\title{
Estimate the pairwise BSDmin for a list of trees
}
\description{
This function calculates the BSDmin distance for a list of trees in parallel. It uses the foreach and doParallel packages. It is a faster version than min.dist.topo.mat when the data set is large. But with small data sets, the parallelization may result in slower computation time.
}
\usage{
min.dist.topo.mat.para(tree.list, para = F, ncore = 1)
}
\arguments{
  \item{tree.list}{
A list containing trees with branch lengths. This can be obtained with the clockstaR function optim.edge.lengths
}
  \item{para}{
Select para = T for parallel computation. In such case select the number of cores with ncore. Default is F for no parallelization
}
  \item{ncore}{
Number of cores to run in parallel mode. Default is 1 for no parallelization
}
}
\details{
None
}
\value{
This function returns a list with two objects. The first is a dist object with the BSDmin between all pairs of trees. The second is a matrix with the scaling factor used for all paris of trees
\item{comp1 }{dist object with minBSD distances among trees}
\item{comp2 }{matrix with scaling factors among trees}
}
\references{
Duchene et al. 
}
\author{
Sebastian Duchene
}
\note{
This function depends on foreach and doParallel packages. Not loading these packages will result in an error. In such case, please use the min.dist.topo.mat function instead.
Please check the user manual and corresponding paper for further details
}

\seealso{
min.dist.topo.mat
Please check the user manual and corresponding paper for further details

}
\examples{
\dontrun{
tr <- rtree(10) # simulate a random tree

trees <- list()
for(i in 1:10){
	trees[[i]] <- tr
	trees[[i]]$edge.length <- abs(rnorm(length(tr$edge.length))) + 0.01 # simulate random branch lengths
}

bsd.matrix.list <- min.dist.topo.mat.para(trees, para = FALSE, ncore = 1)# set para = T and ncore = number of cores to run in parallel

par(mfrow = c(2, 1))
hist(bsd.matrix.list[[1]], main = "BSDmin distances in simulated branch lengths", ylab = "Frequency", xlab = expression(beta))
hist(bsd.matrix.list[[2]], main = "Scaling factors in simulated branch lengths", ylab = "Frequency", xlab = expression(beta))


}


## The function is currently defined as
function (tree.list, para = F, ncore = 1) 
{
    min.dist.topo <- function(tree1, tree2) {
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
            return(min.dist.topo(a, fix.tree))
        })
        return(res)
    }
    if (para == T) {
        cl <- makeCluster(ncore)
        registerDoParallel(cl)
        print(paste("Clusters registered as follows: ", cl))
        res.par <- foreach(s.trees = sub.trees, j = 1:length(tree.list), 
            .packages = c("ape", "phangorn")) \%dopar\% compute.tree.dists(s.trees, 
            tree.list[[j]])
        stopCluster(cl)
    }
    else if (para == F) {
        res.par <- foreach(s.trees = sub.trees, j = 1:length(tree.list), 
            .packages = c("ape", "phangorn")) \%do\% compute.tree.dists(s.trees, 
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
    return(res.list)
  }
}
\keyword{ BSDmin distance }
\keyword{ tree distance score }
