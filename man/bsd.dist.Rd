\name{bsd.dist}
\alias{bsd.dist}
\title{
Estimate sBSDmin distnace between a pair of trees
}
\description{
This function estimates the sBSDmin distance between a pair of trees. This distance metric is described in the original publication. See reference bellow.
}
\usage{
bsd.dist(tree1, tree2)
}
\arguments{
  \item{tree1}{
A phylogenetic tree of class 'phylo'. It should have branch lengths.
}
  \item{tree2}{
A phylogenetic tree of class 'phylo'. It should have branch lengths.

}
}
\details{
See reference
}
\value{
A numeric vector with the sBSDmin (min.bdi.scaled), the scaling factor (scaling.factor), and the unscaled BSD (min.bdi).
}
\references{
Duchêne, S., Molak, M., and Ho, SYW. "ClockstaR: choosing the number of relaxed-clock models in molecular phylogenetic analysis." Bioinformatics (2013): btt665.
}
\author{
Sebastian Duchene
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

# Simulate two trees with identical topologies, but different patterns of among lineage rate variation

tr1 <- rtree(10)

tr2 <- tr1
tr2$edge.length <- sample(tr1$edge.length)

par(mfrow = c(1, 2))
plot(tr1)
plot(tr2)

bsd.dist(tr1, tr2)


## The function is currently defined as
function (tree1, tree2) 
{
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
    tree1.root.scaled$edge.length <- tree1$edge.length * root.scaling
    tree2.root.scaled$edge.length <- tree2.scaled$edge.length * 
        root.scaling
    min.bdi.root.scaled <- dist.topo(tree1.root.scaled, tree2.root.scaled, 
        method = "score")
    res.vect <- c(min.bdi.root.scaled, scaling, min.bdi)
    names(res.vect) <- c("min.bdi.scaled", "scaling.factor", 
        "min.bdi")
    return(res.vect)
  }
}
\keyword{ relaxed-clock }
