\name{cut.trees.beta}
\alias{cut.trees.beta}
\title{
cut.trees.beta calculated the diameter of a dendrogram. If the diameter of the dendrogram is higher than that specified by the user, the dendrogram is cut at its longest edge.
}
\description{
cut.trees.beta calculated the diameter of a dendrogram. If the diameter of the dendrogram is higher than that specified by the user, the dendrogram is cut at its longest edge. This function is internally called by get.all.groups.k and get.all.groups.
}
\usage{
cut.trees.beta(tree, beta = 0.05)
}
\arguments{
  \item{tree}{
This is a dendrogram of the sBSDmin distances among trees. It should be an object of class "phylo"
}
  \item{beta}{
This is the threshold sBSDmin distance. If the diameter of the dendrogram is bellow this value, the dendrogram is cut in two. Otherwise the dendrogram is not cut. The default is 0.05. Note that if the larger than the diameter the dendrogram is not cut.
}
}
\details{
This function is called internally by get.all.groups.k, get.all.groups, and run.clockstar
Please see the user manual for fu}
\value{
\item{comp1 }{A list with the sub dendrograms as objects of class "phylo". If the dendrogram is cut, the lenghth of the list is two, otherwise it will have a length of one with the same dendrogram supplied.}
}
\references{
Duchene et al.
}
\author{
Sebastian Duchene
}
\note{
Please see the user manual for further details
}
\seealso{
get.all.groups.k, get.all.groups, run.clockstar()
}
\examples{
\dontrun{
# generate an example dendrogram
set.seed(1234)
bsd.example <- unroot(rtree(10))

par(mfrow = c(1, 3))
plot(bsd.example, "unrooted")
subdendrograms <- cut.trees.beta(bsd.example, beta = 0.5)
plot(subdendrograms[[1]], "unrooted")
plot(subdendrograms[[2]], "unrooted")
}


## The function is currently defined as
function (tree, beta = 0.05) 
{
    if (length(tree$tip.label) >= 3) {
        tree <- unroot(tree)
    }
    tree.diam <- max(cophenetic(tree))
    pruned.trees <- list()
    if (tree.diam > beta & (length(tree$tip.label) > 2)) {
        di.nodes <- dist.nodes(tree)
        max.edge <- max(tree$edge.length)
        tips <- 1:length(tree$tip.label)
        nodes <- 1:nrow(di.nodes)
        nodes <- nodes[-tips]
        connect.longest.edge <- which(di.nodes == max.edge, arr.ind = T)[1, 
            ]
        if ((connect.longest.edge[1] \%in\% nodes) & (connect.longest.edge[2] \%in\% 
            nodes)) {
            taxa.cut1 <- tips(tree, connect.longest.edge[1])
            taxa.cut2 <- tips(tree, connect.longest.edge[2])
            len.tax <- c(length(taxa.cut1), length(taxa.cut2))
            node.cut <- connect.longest.edge[len.tax == min(len.tax)]
            taxa.cut <- tips(tree, node.cut)
        }
        else if (sum(connect.longest.edge \%in\% nodes) == 1) {
            node.cut <- connect.longest.edge[connect.longest.edge \%in\% 
                nodes]
            taxa.cut <- tree$tip.label[connect.longest.edge[connect.longest.edge \%in\% 
                tips]]
        }
        taxa.prune <- taxa.cut
        if ((length(tree$tip.label) - length(taxa.prune)) >= 
            2) {
            subtree1 <- drop.tip(tree, taxa.prune)
            if (length(subtree1$tip.label) >= 3) {
                subtree1 <- unroot(subtree1)
            }
        }
        else {
            subtree1 <- tree$tip.label[!(tree$tip.label \%in\% 
                taxa.prune)]
        }
        if (length(taxa.prune) >= 2) {
            subtree2 <- drop.tip(tree, tree$tip.label[!(tree$tip.label \%in\% 
                taxa.prune)])
            if (length(subtree2$tip.label) >= 3) {
                subtree2 <- unroot(subtree2)
            }
        }
        else {
            subtree2 <- taxa.prune
        }
        pruned.trees[[1]] <- subtree1
        pruned.trees[[2]] <- subtree2
    }
    else if (tree.diam < beta) {
        pruned.trees[[1]] <- tree
    }
    else {
        pruned.trees[[1]] <- tree$tip.label[1]
        pruned.trees[[2]] <- tree$tip.label[2]
    }
    return(pruned.trees)
  }
}

\keyword{ sBSDmin }
\keyword{ diameter of a graph }
