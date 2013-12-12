\name{min.dist.topo}
\alias{min.dist.topo}
\title{
This function estimates the sBSDmin value described in Duchene et al. 
}
\description{
The function first obtains the minimum branch distance score between a pair of trees by using univariate optimisation. The trees are then rescaled so that the mean sum of their length is = 0.05. 
}
\usage{
min.dist.topo(tree1, tree2)
}
\arguments{
  \item{tree1}{
This is a phylogenetic tree as an object of class 'phylo' with branch lengths.
}
  \item{tree2}{
This object is a phylogenetic tree as an object of class phylo with branch lengths.
}
}
\details{
In the method described in Duchene et al. the topologies of tree1 and tree2 are identical, and they differ in their branch lengths.
}
\value{
The function returns an object of class 'numeric' with three elements as follows: min.bdi.scaled is the beta value for the pair of trees
scaling.factor is the scaling factor so that the mean sum of the tree lengths is 0.05.
min.bdi is the BSDmin score proposed by Kuhner and Felsenstein (1994) with modifications by Duchene et al.  
}
\references{
See Duchene et al. and Kuhner and Felsenstein (1994)}
\author{
Sebastian Duchene
}
\note{
Please see the corresponding user manual and paper for details of the implementation (Duchene et al. )
}
\seealso{
min.dist.topo.mat and min.dist.topo.mat.para estimate the pairwise beta for a list of trees.
}
\examples{
\dontrun{
# Simulate a tree and store it in two objects
tr1 <- rtree(50)
tr2 <- tr1
# change the branch lengths of one of the trees 
tr2$edge.length <- abs(rnorm(length(tr1$edge.length), 0.02, 2)) 

# plot the trees to verify that their specific are different 
par(mfrow = c(1, 2))
plot(tr1, show.tip.label = FALSE, direction = "right")
plot(tr2, show.tip.label = FALSE, direction = "left")

#Estimate the beta
min.dist.topo(tr1, tr2)
}

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
