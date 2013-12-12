\name{diagnostics.clockstar}
\alias{diagnostics.clockstar}
\title{
diagnostics.clockstar
}
\description{
diagnostics.clockstar is a function to inspect the results from the ClockstaR analysis.
}
\usage{
diagnostics.clockstar(groups.obj, save.plots = FALSE, plots.file = "clockstar.diagnostics.pdf", interactive = TRUE)
}
\arguments{
  \item{groups.obj}{
This should be the output from calling function get.all.groups. For instance, use:
example.data <- get.all.groups(sbsdmin.matrix, pam.results = T)
example.diagnostics <- diagnostics.clockstar(example.data)

sbsdmin.matrix is the matrix of the pairwise sBSDmin matrix estimated in ClockstaR.
}
  \item{save.plots}{
Should the diagnostic plots should be saved to a pdf file? If so, select T.
}
  \item{plots.file}{
The name for the pdf file with the plots.
}
  \item{interactive}{
If this argument is TRUE then ClockstaR will print the diagnostic plots in the active graphics device.
}

}
\details{
This function can be used to assess the partitioning strategies for different values of k. Please see the tutorial for an example.
}
\value{
\item{comp1 }{This function prints the Gap statistic for different values of k, and an image matrix with colours to represent the partition assignment. 
}
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
# Simulate 10 data sets with three different clocks
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

bsd.matrix.list <- as.matrix(min.dist.topo.mat(tree.list)[[1]])

example.partitions <- get.all.groups(bsd.matrix.list, pam.results = TRUE)

example.diagnostics <- diagnostics.clockstar(example.partitions)
}
## The function is currently defined as
function(groups.obj, save.plots = FALSE, plots.file = "clockstar.diagnostics.pdf", interactive = TRUE)
{
    if(length(groups.obj) == 3){
        stop("The object specified is not the correct format. Please provide an the output of function get.all.groups. Please use pam.restuls = T")
    }
    gap.stat <- groups.obj[[2]][[1]][, 3]
    ks <- 1:length(gap.stat)
    errs.up <- gap.stat + groups.obj[[2]][[1]][, 4]
    errs.down <- gap.stat - groups.obj[[2]][[1]][, 4]
    k.matrix <- matrix(NA, nrow = nrow(groups.obj[[4]]), ncol = 1)
    for(i in 1:length(ks)){
        part <- pam(groups.obj[[4]], k = ks[i])
        k.matrix <- cbind(k.matrix, as.matrix(part[[3]]))
    }
    k.matrix <- k.matrix[, 2:ncol(k.matrix)]
    plot.matrices <- function(clust.matrix){
        par(mar = c(2, 4, 2, 4))
        image(t(as.matrix(clust.matrix)), axes = F, col = sample(rainbow(ncol(clust.matrix)*100), ncol(clust.matrix)), ylab = "Data subset")
        mtext(text = rownames(clust.matrix), side = 2, line = 0.3, at = seq(0, 1, 1/(nrow(clust.matrix) - 1)), las = 1, cex = 0.6)
        mtext(text = paste0("k=", ks), side = 3, line = 0.3, at = seq(0, 1, 1/(ncol(clust.matrix) - 1)), las = 1, cex = 0.6)
    }
    if(interactive == TRUE){
    	par(mfrow = c(2, 1))
    	par(mar = c(4, 4, 4, 4))
    	errbar(x = ks, y = gap.stat, yplus = errs.up, yminus = errs.down, type = "b", pch = 20, xlab = expression(italic(k)), ylab = "Gap statistic")
    	points(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
    	lines(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
    	lines(1:(length(groups.obj[[3]]) + 1), y = rep(0, length(groups.obj[[3]]) + 1), col = "blue")
    	legend(x = max(ks)*0.7, y = min(gap.stat)*0.5, legend = c("Gap statistic and SE" ,"Delta Gap statistic", "Gap statistic = 0"), cex = 0.8, fill = c("black", "red", "blue"), bty = "n")
    	plot.matrices(k.matrix)
    }
    if(save.plots == TRUE){
        pdf(file = plots.file)
    	par(mfrow = c(2, 1))
        par(mar = c(4, 4, 4, 4))
        errbar(x = ks, y = gap.stat, yplus = errs.up, yminus = errs.down, type = "b", pch = 20, xlab = expression(italic(k)), ylab = "Gap statistic")
        points(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
        lines(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
        lines(1:(length(groups.obj[[3]]) + 1), y = rep(0, length(groups.obj[[3]]) + 1), , col = "blue")
        legend(x = max(ks)*0.7, y = min(gap.stat)*0.5, legend = c("Gap statistic and SE" ,"Delta Gap statistic", "Gap statistic = 0"), cex = 0.8, fill = c("black", "red", "blue"), bty = "n")
        plot.matrices(k.matrix)
        dev.off()
    }
    return(k.matrix)
}
}

\keyword{ sBSDmin }
\keyword{ diameter of a graph }
