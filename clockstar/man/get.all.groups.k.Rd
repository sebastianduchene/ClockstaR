\name{get.all.groups.k}
\alias{get.all.groups.k}
\title{
This function finds the best partitioning strategy for a fixed number of partitions, as described in Duchene et al.
}
\description{
This function uses the PAM algorithm to assign data subsets to a fixed number of partitions (k). The input is a matrix or a data frame with the sBSDmin distance between all pairs of trees in the data set.
If an object of dendrogram of sBSDmin distances (class phylo) is specified, instead of a matrix or data frame, the function cuts the dendrogram at its longest edge, producing two dendrograms. The process is repeated until a k number of groups is obtained. This function is used when clockstar is run in k mode.
}
\usage{
get.all.groups.k(data.obj, k.man = 1, save.partitions = F, file.name = "partitions.txt")
}
\arguments{
  \item{data.obj}{
This is the data object. It can be the sBSDmin distance among all pairs of trees, or a dendrogram constructed using these distances.

If the data object is the sBSDmin distances, then it should correspond to a data.frame or a matrix object. In this case the partitions are defined using the PAM algorithm as implemented in the cluster package.

If the data object is a dendrogram, the function cuts it along its longest edge, producing two dendrograms. The process is repeated until the diameter of all dendrograms is the selected k, as specified in the parameter k.man.
}
  \item{k.man}{
The selected number of partitions
}
  \item{save.partitions}{
To save the partitions in a text file select T. The default is F
}
  \item{file.name}{
This is the name for the output file with the partitions. This parameter is only used when save.partitions = T.
}
}
\details{
See Duchene et al. and the user manual for more details
}
\value{
A list where each item is the data subsets in the partition. Partitions are named "partition_1" through "partition_k"
}
\references{
Duchene et al. 
}
\author{
Sebastian Duchene
}

\seealso{
get.all.groups is similar to get.all.groups.k but in get.all.groups the number of partitions is selected automatically using the GAP statistic and the PAM algorithm as implemented in the cluster packages.
See user manual and Duchene et al. for futher details
}
\examples{
# This is the same example in the package description
# generate a list of trees with three patterns of among-lineage rate variation
\dontrun{
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

groups.example <- get.all.groups.k(sbsdmin, k.man = 3)

#print the partitions

print(groups.example)
	
}
## The function is currently defined as
function (data.obj, k.man = 1, save.partitions = F, file.name = "partitions.txt") 
{
    if (class(data.obj) \%in\% c("data.frame", "matrix")) {
        print("The data object is a data frame or a matrix. The partitions will be defined with the PAM algorithm")
        if (k.man == 1) {
            stop("The selected number of partitions is 1. Please select k.man for  >1")
        }
        else {
            pam.k <- pam(data.obj, k = k.man)
            partitions <- unique(pam.k$clustering)
            res.data <- list()
            for (i in 1:length(partitions)) {
                res.data[[i]] <- names(which(pam.k$clustering == 
                  partitions[i]))
            }
            names(res.data) <- paste("Partition_", 1:length(res.data))
            if (save.partitions == T) {
                cat(paste("Partitions selected with manual mode for k =", 
                  k.man, "\n"), file = file.name)
                for (m in 1:length(res.data)) {
                  cat(names(res.data[m]), file = file.name, sep = "\n", 
                    append = T)
                  cat(res.data[[m]], file = file.name, append = T)
                  cat("\n", file = file.name, append = T)
                }
            }
            return(res.data)
        }
    }
    else if (class(data.obj) == "phylo") {
        print("The data object is a dendrogram (phylo object). The partitions will be defined but cutting the dendrobran along the longest edge")
        tree <- data.obj
        k = 2
        tree.list <- list()
        tree.list[[1]] <- tree
        while (length(tree.list) < k) {
            diams <- sapply(tree.list, function(tr) {
                if (class(tr) == "phylo") {
                  return(max(cophenetic(tr)))
                }
                else {
                  return(0)
                }
            })
            tree.to.cut <- tree.list[[which(diams == max(diams))]]
            tree.list <- tree.list[-which(diams == max(diams))]
            tree.list[c(length(tree.list) + 1, length(tree.list) + 
                2)] <- cut.trees.beta(tree.to.cut, beta = 0)
        }
        for (l in 1:length(tree.list)) {
            if (class(tree.list[[l]]) == "phylo") {
                tips <- tree.list[[l]]$tip.label
                tree.list[[l]] <- tips
            }
        }
        names(tree.list) <- paste("Partition_", 1:length(tree.list))
        if (save.partitions == T) {
            cat(paste("partitions with selected k =", k, "\n"), 
                file = file.name)
            for (m in 1:length(tree.list)) {
                cat(names(tree.list[m]), file = file.name, sep = "\n", 
                  append = T)
                cat(tree.list[[m]], file = file.name, append = T)
                cat("\n", file = file.name, append = T)
            }
        }
        return(tree.list)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
