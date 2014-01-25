\name{get.all.groups}
\alias{get.all.groups}
\title{
This function finds the best number of partitions (k) for a data set, and the data subsets in each partitions, as described in Duchene et al.
}
\description{
This function uses the GAP statistic with the PAM algorithm. The input is a matrix or a data frame with the sBSDmin distance between all pairs of trees in the data set.
}
\usage{
get.all.groups(data.obj, n.b = 500, pam.results = F, save.partitions = F, file.name = "partitions.txt", beta = 0.03, ...)
}
\arguments{
  \item{data.obj}{
This is the data object. It can be the sBSDmin distance among all pairs of trees, or a dendrogram constructed using these distances.

If the data object is the sBSDmin distances, then it should correspond to a data.frame or a matrix object. In this case the partitions are defined using the GAP statistic with the PAM algorithm as implemented in the cluster package.

If the data object is a dendrogram, the function cuts it along its longest edge, producing two dendrograms. The process is repeated until the diameter of all dendrograms is < a threshold value (known as beta). Note that this method is deprecated as the PAM as the choice of beta is subjective. Only use for test purposes.
}
  \item{n.b}{
Number of bootstrap replicates to estimate the GAP statistic.
}
  \item{pam.results}{
Select T to obtain the results from the GAP statistics. This would be a data frame as estimated in the cluster package. For more details use ?clusGap. If pam.results=T, the function returns a list where the first element is a list of the partitions, and the second is the data frame with the results from clusGap. The default is F.
}
  \item{save.partitions}{
To save the partitions in a text file select T. The default is F
}
  \item{file.name}{
This is the name for the output file with the partitions. This parameter is only used when save.partitions = T.
}
  \item{beta}{
In most analyses this parameter will not be used. This is the threshold value for cutting the dendrogram successively, as described in the above section for data.obj. Please note that this is deprecated, and should only be used for test purposes.
}
  \item{\dots}{
Additional arguments for clusGap.
}
}
\details{
See Duchene et al. and the tutorial for more details.
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
\note{
get.all.groups.k is similar to get.all.groups but the user can select the number of partitions (k) instead of the the automatic selection with the GAP statistic implemented here.
See user manual and Duchene et al. for futher details
}

\seealso{
Use help(get.all.groups.k) for a version of the function where the number of partitions can be selected.
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

groups.example <- get.all.groups(sbsdmin)

#print the partitions

print(groups.example)
}	


## The function is currently defined as
function (data.obj, n.b = 500, pam.results = F, save.partitions = F, 
    file.name = "partitions.txt", beta = 0.03, ...) 
{
    if (class(data.obj) \%in\% c("data.frame", "matrix")) {
        kmax <- nrow(data.obj) - 1
        print("FINDING THE BEST NUMBER OF PARTITIONS (k) WITH GAP STATISTIC AND THE PAM ALGORITHM")
        print(kmax)
        pam.gap <- clusGap(data.obj, pam, B = n.b, K.max = kmax)
        rownames(pam.gap[[1]]) <- paste0("k=", 1:nrow(pam.gap[[1]]))
        gaps <- pam.gap$Tab[, 3]
        diffs <- vector()
        for (i in 1:(length(gaps) - 1)) {
            diffs <- c(diffs, gaps[i] - gaps[i + 1])
        }
        names(diffs) <- paste0("k=", 1:length(diffs))
        all.k <- (1:length(diffs))[diffs > 0]
        best.k <- maxSE(pam.gap[[1]][, 3], pam.gap[[1]][, 4])
        pam.k <- pam(data.obj, k = best.k, ...)
        partitions <- unique(pam.k$clustering)
        res.data <- list()
        for (i in 1:length(partitions)) {
            res.data[[i]] <- names(which(pam.k$clustering == 
                partitions[i]))
        }
        names(res.data) <- paste("Partition_", 1:length(res.data))
        if (save.partitions == T) {
            cat("Partitions selected with automatic mode \n", 
                file = file.name)
            for (m in 1:length(res.data)) {
                cat(names(res.data[m]), file = file.name, sep = "\n", 
                  append = T)
                cat(res.data[[m]], file = file.name, append = T)
                cat("\n", file = file.name, append = T)
            }
        }
        if (pam.results == T) {
            lis.res <- list(res.data, pam.gap, diffs, data.obj)
            return(lis.res)
        }
        else if (pam.results == F) {
            return(res.data)
        }
    }
    else if (class(data.obj) == "phylo") {
        tree <- data.obj
        temp.list <- list()
        min.list <- list()
        temp.list[[1]] <- tree
        get.diameter <- function(tr) {
            if (class(tr) == "phylo") {
                return(max(cophenetic(tr)))
            }
            else {
                return(0)
            }
        }
        diams <- sapply(temp.list, get.diameter)
        if (any(diams > beta)) {
            while (length(temp.list) != 0) {
                diams <- sapply(temp.list, get.diameter)
                if (sum(diams <= beta) > 0) {
                  diams.beta <- seq(from = 1, to = length(diams))[diams <= 
                    beta]
                  for (i in diams.beta) {
                    min.list[[length(min.list) + 1]] <- temp.list[[i]]
                  }
                }
                else {
                  diams.beta = 0
                }
                if (sum(diams > beta) > 0) {
                  diams.non.beta <- seq(from = 1, to = length(diams))[diams > 
                    beta]
                  temp.list.non.beta <- list()
                  for (j in 1:length(diams.non.beta)) {
                    temp.list.non.beta[[j]] <- temp.list[[diams.non.beta[j]]]
                  }
                  temp.list <- temp.list.non.beta
                  sub.list <- list()
                  for (k in 1:length(temp.list)) {
                    cut.temp <- cut.trees.beta(temp.list[[k]], 
                      beta)
                    sub.list[[length(sub.list) + 1]] <- cut.temp[[1]]
                    if (length(cut.temp) == 2) {
                      sub.list[[length(sub.list) + 1]] <- cut.temp[[2]]
                    }
                  }
                  temp.list = sub.list[1:length(sub.list)]
                }
                else {
                  temp.list <- list()
                }
            }
        }
        else {
            min.list <- temp.list
        }
        for (l in 1:length(min.list)) {
            if (class(min.list[[l]]) == "phylo") {
                tips <- min.list[[l]]$tip.label
                min.list[[l]] <- tips
            }
        }
        names(min.list) <- paste("Partition_", 1:length(min.list))
        if (save.partitions == T) {
            cat(paste("partitions with selected beta =", beta, 
                "\n"), file = file.name)
            for (m in 1:length(min.list)) {
                cat(names(min.list[m]), file = file.name, sep = "\n", 
                  append = T)
                cat(min.list[[m]], file = file.name, append = T)
                cat("\n", file = file.name, append = T)
            }
        }
        return(min.list)
    }
    else {
        stop("Please supply the distance between trees. This should be an object of classs Matrix of data.frame. For objects of class dist, use as.matrix(data) to change the class of the object")
    }
  }
}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
