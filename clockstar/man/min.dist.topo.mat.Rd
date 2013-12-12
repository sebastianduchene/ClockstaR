\name{min.dist.topo.mat}
\alias{min.dist.topo.mat}
\title{
Estimate the pairwise BSDmin for a list of trees
}
\description{
This function calculates the BSDmin distance for a list of trees. For a parallel version use min.dist.topo.mat.para
}
\usage{
min.dist.topo.mat(tree.list)
}
\arguments{
  \item{tree.list}{
A list containing trees with branch lengths. This can be obtained with the clockstaR function optim.edge.lengths
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
Please see the user manual and corresponding paper for further details
}

\seealso{
min.dist.topo.mat.para
Please check the user manual and corresponding paper for further details
}
\examples{
\dontrun{
trees <- list()
tr <- rtree(10)
for(i in 1:10){
	trees[[i]] <- tr
	trees[[i]]$edge.length <- abs(rnorm(length(tr$edge.length))) + 0.01 # simulate random branch lengths
}

bsd.matrix.list <- min.dist.topo.mat(trees)

par(mfrow = c(2, 1))
hist(bsd.matrix.list[[1]], main = "BSDmin distances in simulated branch lengths", ylab = "Frequency", xlab = expression(beta))
hist(bsd.matrix.list[[2]], main = "Scaling factors in simulated branch lengths", ylab = "Frequency", xlab = expression(beta))
}

## The function is currently defined as
function (tree.list) 
{
    d.mat <- matrix(NA, nrow = length(tree.list), ncol = length(tree.list))
    rownames(d.mat) <- names(tree.list)
    colnames(d.mat) <- names(tree.list)
    s.mat <- d.mat
    print("Estimating tree distances")
    if (length(tree.list) > 2) {
        d.mat.lin <- vector()
        d.mat.lin <- sapply(2:nrow(d.mat), function(a) {
            print(paste("estimating distances", a - 1, "of", 
                nrow(d.mat) - 1))
            lapply(tree.list[1:(a - 1)], function(y) {
                min.dist.topo(tree1 = y, tree2 = tree.list[[a]])
            })
        })
        for (a in 1:length(d.mat.lin)) {
            vec.temp.dist <- vector()
            vec.temp.scale <- vector()
            for (b in 1:length(d.mat.lin[[a]])) {
                vec.temp.dist[b] <- d.mat.lin[[a]][[b]][1]
                vec.temp.scale[b] <- d.mat.lin[[a]][[b]][2]
            }
            d.mat[a + 1, 1:length(vec.temp.dist)] <- vec.temp.dist
            s.mat[a + 1, 1:length(vec.temp.dist)] <- vec.temp.scale
        }
    }
    else if (length(tree.list) == 2) {
        d.mat.lin <- min.dist.topo(tree.list[[1]], tree.list[[2]])
        d.mat[2, 1] <- d.mat.lin[1]
        s.mat[2, 1] <- d.mat.lin[2]
    }
    else {
        print("there was only one tree in the list or the objects are not of class 'phylo'")
    }
    res.list <- list()
    res.list[[1]] <- as.dist(d.mat)
    res.list[[2]] <- s.mat
    return(res.list)
  }
}
\keyword{ BSDmin distance }
\keyword{ tree distance score }
