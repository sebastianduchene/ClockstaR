\name{clockstar-package}
\alias{clockstar-package}
\alias{clockstar}
\docType{package}
\title{
ClockstaR: Selecting partitioning strategies for relaxed phylogenetic analysis
}
\description{
This package implements the method described in Duchene et al. It assigns multigene data subsets to molecular clock partitions in relaxed molecular clock analysis.
}
\details{
\tabular{ll}{
Package: \tab clockstar\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-02\cr
License: \tab What license is it under?\cr
}
Please see the user manual for a tutorial on how to run ClockstaR
}
\author{
Sebastian Duchene

Maintainer: Sebastian Duchene <sebastian.duchene@sydney.edu.au>
}
\references{
Duchene et al. 
Duchene et al.
}
\keyword{  Relaxed molecular clock  }
\seealso{
Please visit the website bellow for more tutorials and more information
http://sydney.edu.au/science/biology/meep/software/
}
\examples{
# This is an example of using ClockstaR in automatic mode.

# generate a list of trees with three patterns of among-lineage rate variation

fixed.topology <- rtree(10)
fixed.topology$edge.length <- NULL
tree.list <- list()
set.seed(1234)
patterns1 <- abs(rnorm(18, 0, 10))
set.seed(3456)
patterns2 <- abs(rnorm(18, 0, 10))
set.seed(0987)
patterns3 <- abs(rnorm(18, 0, 10))
for(i in 1:3){
	tree.list[[i]] <- fixed.topology
	tree.list[[i]]$edge.length <- abs(patterns1 + rnorm(18))
}

for(i in 4:7){
	tree.list[[i]] <- fixed.topology
	tree.list[[i]]$edge.length <- abs(patterns2 + rnorm(18))
}

for(i in 8:10){
	tree.list[[i]] <- fixed.topology
	tree.list[[i]]$edge.length <- abs(patterns3 + rnorm(18))
}

# to inspect some of the trees:
par(mfrow = c(3, 1))
plot(tree.list[[1]])
plot(tree.list[[4]])
plot(tree.list[[10]])

# estimate sBSDmin for all pairs of trees

sbsdmin <- as.matrix(min.dist.topo.mat(tree.list)[[1]])

groups.example <- get.all.groups(sbsdmin)

#print the partitions

print(groups.example)
	
	

}
