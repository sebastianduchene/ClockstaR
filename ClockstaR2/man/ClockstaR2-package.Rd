\name{ClockstaR2-package}
\alias{ClockstaR2-package}
\alias{ClockstaR2}
\docType{package}
\title{
ClockstaR: choosing the number of relaxed-clock models in molecular phylogenetic analysis
}
\description{
Relaxed molecular clocks allow the phylogenetic estimation of evolutionary timescales even when substitution rates vary among branches. In analyses of large multigene datasets, it is often appropriate to use multiple relaxed-clock models to accommodate differing patterns of rate variation among genes. We present ClockstaR, a method for selecting the number of relaxed clocks for multigene datasets.
}
\details{
\tabular{ll}{
Package: \tab ClockstaR2\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-03-21\cr
License: \tab GNU\cr
}
}
\author{
Sebastian Duchene
Maintainer: Sebastian Duchene <sebastian.duchene@sydney.edu.au>
}
\references{
Duchêne, Sebastián, Martyna Molak, and Simon YW Ho. "ClockstaR: choosing the number of relaxed-clock models in molecular phylogenetic analysis." Bioinformatics (2013): btt665.
}

\keyword{ molecular-clock }
\seealso{
NONE.
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

# Estimate sBSDmin distance for all pairs of trees:
trees.bsd <- bsd.matrix(trees.list)

# Find the optimal number of partitions:
partitions.object <- partitions(trees.bsd)

# plot partitions in two graphics devices
plot(partitions.object)


}
