\name{partitions}
\alias{partitions}
\title{
function partitions applied to objects of class bsd
}
\description{
This is the definition of method partitions for bsd objects. See partitions.bsd for more details.
}
\usage{
partitions(bsd.object, ...)
}
\arguments{
\item{bsd.object}{An object of class bsd, obtained with bsd.matrix or bsd.matrix.para}
}
\author{
Sebastian Duchene
}

\seealso{
partitions.bsd
}
\examples{

\dontrun{

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

partitions(trees.bsd)


}

## The function is currently defined as
function (bsd.object, ...) 
UseMethod("partitions")
}
