\name{bsd.matrix}
\alias{bsd.matrix}
\title{
Estimate the sBSDmin distance for all pairs of trees in a list
}
\description{
bsd.matrix estimates the sBSDmin distance for all paris of trees in a list. It returns an object of class bsd with the pairswise distances and the scaling factors. The scaling factors represent the differences in the magnitude of the rates among trees. See the help file for bsd.dist for more details.
}
\usage{
bsd.matrix(tree.list)
}
\arguments{
  \item{tree.list}{
A list of trees. It can be a list where each item is a tree, or an object of class multiPhylo. The trees can have names. See the example bellow.
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


## The function is currently defined as
function (tree.list) 
{
    d.mat <- matrix(NA, nrow = length(tree.list), ncol = length(tree.list))
    rownames(d.mat) <- names(tree.list)
    colnames(d.mat) <- names(tree.list)
    s.mat <- d.mat
    print("Estimating tree distances")
    if (length(tree.list) > 3) {
        d.mat.lin <- vector()
        d.mat.lin <- sapply(2:nrow(d.mat), function(a) {
            print(paste("estimating distances", a - 1, "of", 
                nrow(d.mat) - 1))
            lapply(tree.list[1:(a - 1)], function(y) {
                bsd.dist(tree1 = y, tree2 = tree.list[[a]])
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
    else if (length(tree.list) <= 3) {
        stop("The number of gene trees is <= 3. ClockstaR requires at least gene 4 trees")
    }
    d.mat.lin <- bsd.dist(tree.list[[1]], tree.list[[2]])
    d.mat[2, 1] <- d.mat.lin[1]
    s.mat[2, 1] <- d.mat.lin[2]
    res.list <- list()
    res.list[[1]] <- as.dist(d.mat)
    res.list[[2]] <- s.mat
    class(res.list) <- "bsd"
    return(res.list)
  }
}
\keyword{ sBSDmin }
