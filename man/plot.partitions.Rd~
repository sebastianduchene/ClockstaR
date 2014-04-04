\name{plot.partitions}
\alias{plot.partitions}
\title{
plot objects of class partitions
}
\description{
This function plots the results from the gap statistic from package cluster, and a matrix with colours for the partition assignments. 
}
\usage{
plot.partitions(partitions.object, save.plot = F, file.name = "results_clockstar")
}
\arguments{
  \item{partitions.object}{
An object of class 'partitions', obtained with partitions.bsd
}
  \item{save.plot}{
logical. Select T to save the plots to a pdf file, or F to print in two new graphics devices.
}
  \item{file.name}{
A character. If save.plot == T, the name of the file to print the results.
}
}
\author{
Sebastian Duchene
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






## The function is currently defined as
function (partitions.object, save.plot = F, file.name = "results_clockstar") 
{
    image(t(partitions.object$parts.mat), col = sample(rainbow(ncol(partitions.object$parts.mat))), 
        axes = F, ylab = "Data subset", main = "Partitioning scheeme for values of k")
    mtext(text = rownames(partitions.object$parts.mat), side = 2, 
        line = 0.3, at = seq(0, 1, 1/(nrow(partitions.object$parts.mat) - 
            1)), las = 1, cex = 0.6)
    mtext(text = colnames(partitions.object$parts.mat), side = 3, 
        line = 0.3, at = seq(0, 1, 1/(ncol(partitions.object$parts.mat) - 
            1)), las = 1, cex = 0.6)
    if (save.plot) {
        dev.copy2pdf(file = paste0(file.name, "_matrix.pdf"))
    }
    dev.new()
    plot(partitions.object$clusterdata)
    if (save.plot) {
        dev.copy2pdf(file = paste0(file.name, "_gapstats.pdf"))
        sapply(dev.list(), function(x) dev.off(x))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
