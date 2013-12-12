\name{optim.edge.lengths}
\alias{optim.edge.lengths}
\title{
This function estimates the branch lengths of a fixed tree topology for multiple alignment data sets stored in a folder.
}
\description{
The function uses the Maximum Likelihood approach implemented in the phangorn package to optimise the branch lengths of a tree for different data sets. When model.test = T the best model is selected according to the BIC. 
}
\usage{
optim.edge.lengths(directory, fixed.tree, form = "fasta", model.test = F, save.trees = F, tree.file.names = "output", para = F, ncore = 1)
}
\arguments{
  \item{directory}{
This is the directory with the individual alignments in fasta or nexus format.
}
  \item{fixed.tree}{
This is the tree for which the branch lengths are estimated. In this function the tree should be an R object. In the run.clockstar version, the tree can be stored in newick format in the folder with the alignments and does not need to be loaded previously in R.
}
  \item{form}{
The format of the alignments. It can be "nexus", "fasta", or any abbreviation of these. The default is "fasta"
}
  \item{model.test}{
When mode.test = T, the best fitting substitution model is obtained for each data set using the BIC. Otherwise the branch length optimisation is performed with the Jukes-Cantor model.
}
  \item{save.trees}{
If the individual trees should be saved, this option should be set to save.trees = T. In this case the trees are written in a newick format. Note that the topologies are identical, but the branch lengths are those of each data set. The trees are saved in the order that the data sets are read into R, which can be verified in the the second item of the list returned by the function. The default is save.trees = F
}
  \item{tree.file.names}{
The file name for the optimised trees.
}
  \item{para}{
In multicore machines with the foreach and doParallel packages the function can be run in parallel. Use para = T to enable parallelisation. The default is para = F.
}
  \item{ncore}{
Number of cores to run in parallel. Only meaningful if para = T
}
}
\details{
This function is called internally when run.clockstar() is used.
}
\value{
\item{comp1 }{The first item of the list is a list of trees with the branch lengths optimised for each alignment. Although this is a list, it is defined with class multiPhylo for ease of use with other phylogenetic R packages. The names of the individual trees are those of the data sets.
}
\item{comp1 }{The second item of the list is a matrix with the models selected for each data set and their corresponding BIC score. When optim.edge.lengths is used internally from the run.clockstar() function the matrix is saved in the results folder. If para = T, the matrix will be filled with NA, even if model.test = T. Although the parellel version performs model testing, the models selected according to the BIC are not saved. 
}
}
\references{
Duchene et al.
}
\author{
Sebastian Duchene
}
\note{
See the user manual for ClockstaR
}
\seealso{
run.clockstar, which uses this function internally for a more user friendly implementation.
}
\examples{
\dontrun{
system("mkdir optim_edge_lengths_test")
tr <- rtree(10)


setwd("./optim_edge_lengths_test")
for(i in 1:10){
	tr.temp <- tr
	tr.temp$edge.length <- runif(length(tr$edge.length))
	data.temp <- as.DNAbin(simSeq(tr))
	write.dna(data.temp, file = paste("sim_data", i, ".fasta", sep = ""), format = "fasta")
}
write.tree(tr, file = "sim_fix_tree.tre")

setwd("..")

optim.data <- optim.edge.lengths(directory = "optim_edge_lengths_test", fixed.tree = tr, model.test = TRUE)

par(mfrow = c(2, 5))
par(mar = c(4, 4, 4, 4))
for(i in 1:10) hist(optim.data[[1]][[i]]$edge.length, col = i, main = "", ylab = "Frequency", xlab = paste("Branch lengths for", names(optim.data[[1]])[i]))

dev.new()
par(mfrow = c(2, 5))
par(mar = c(4, 4, 4, 4))
for(i in 1:10){
	plot(optim.data[[1]][[i]])
	axisPhylo()
}

print(optim.data[[2]])

}
## The function is currently defined as
function (directory, fixed.tree, form = "fasta", model.test = F, 
    save.trees = F, tree.file.names = "output", para = F, ncore = 1) 
{
    options(warn = -1)
    directory = paste(directory, "/", sep = "")
    file.names <- dir(directory)
    drop.tip <- ape::drop.tip
    file.names <- file.names[grepl(form, file.names)]
    model.table <- matrix(NA, nrow = length(file.names), ncol = 3)
    colnames(model.table) <- c("file", "BIC", "model")
    model.table[, 1] <- file.names
    data.files <- list()
    trees.opt <- list()
    print("reading files")
    for (a in 1:length(file.names)) {
        data.files[[a]] <- read.dna(paste(directory, file.names[a], 
            sep = ""), format = form)
    }
    if (para == F) {
        for (b in 1:length(file.names)) {
            tax.keep.temp <- fixed.tree$tip.label \%in\% rownames(data.files[[b]])
            trees.opt[[b]] <- drop.tip(fixed.tree, as.character(fixed.tree$tip.label[!tax.keep.temp]))
            trees.opt[[b]]$edge.length <- rtree(nrow(data.files[[b]]))$edge.length
            pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]), 
                inv = 0, shape = 1, k = 1)
            print(paste("model testing dataset", file.names[b], 
                b, "of", length(file.names)))
            if (model.test == T) {
            	dat.temp <- phyDat(data.files[[b]])
				model.temp <- modelTest(dat.temp, multicore = T)
                model.table[b, 2:3] <- c(model.temp$BIC[model.temp$BIC == 
                  min(model.temp$BIC)], model.temp$Model[model.temp$BIC == 
                  min(model.temp$BIC)])
                best.model.temp <- model.temp$Model[model.temp$BIC == 
                  min(model.temp$BIC)][1]
            }
            else if (model.test == F) {
                best.model.temp <- "JC"
            }
            if (length(grep("+G", best.model.temp)) == 0 && length(grep("+I", 
                best.model.temp)) == 0) {
                pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]))
                trees.opt[[b]] <- optim.pml(pml.temp, optEdge = T)$tree
            }
            else if (length(grep("+G", best.model.temp)) == 1 && 
                length(grep("+I", best.model.temp)) == 0) {
                pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]), 
                  optInv = T)
                trees.opt[[b]] <- optim.pml(pml.temp, optEdge = T, 
                  optGamma = T)$tree
            }
            else if (length(grep("+G", best.model.temp)) == 0 && 
                length(grep("+I", best.model.temp)) == 1) {
                pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]), 
                  optGamma = T)
                trees.opt[[b]] <- optim.pml(pml.temp, optEdge = T, 
                  optInv = T)$tree
            }
            else if (length(grep("+G", best.model.temp)) == 1 && 
                length(grep("+I", best.model.temp)) == 1) {
                pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]), 
                  optInv = T, optGamma = T)
                trees.opt[[b]] <- optim.pml(pml.temp, optEdge = T, 
                  optGamma = T, optInv = T)$tree
            }
            print(paste("optimized edge lengths for tree", b, 
                "of", length(file.names)))
        }
    }
    optim.trees.par <- function(dat) {
        tree.par <- fixed.tree
        dat.file <- dat
        tax.keep.temp <- fixed.tree$tip.label \%in\% rownames(dat.file)
        tree.par <- drop.tip(fixed.tree, as.character(fixed.tree$tip.label[!tax.keep.temp]))
        lens.temp <- rtree(nrow(dat.file))
        tree.par$edge.length <- lens.temp$edge.length
        pml.temp <- pml(tree.par, phyDat(dat.file), inv = 0, 
            shape = 1, k = 1)
        if (model.test == T) {
        	dat.temp <- phyDat(dat)
			model.temp <- modelTest(dat.temp, multicore = T)
            best.model.temp <- model.temp$Model[model.temp$BIC == 
                min(model.temp$BIC)][1]
        }
        else if (model.test == F) {
            best.model.temp <- "JC"
        }
        if (length(grep("+G", best.model.temp)) == 0 && length(grep("+I", 
            best.model.temp)) == 0) {
            pml.temp <- pml(tree.par, phyDat(dat.file))
            tree.par <- optim.pml(pml.temp, optEdge = T)$tree
        }
        else if (length(grep("+G", best.model.temp)) == 1 && 
            length(grep("+I", best.model.temp)) == 0) {
            pml.temp <- pml(tree.par, phyDat(dat.file), optInv = T)
            tree.par <- optim.pml(pml.temp, optEdge = T, optGamma = T)$tree
        }
        else if (length(grep("+G", best.model.temp)) == 0 && 
            length(grep("+I", best.model.temp)) == 1) {
            pml.temp <- pml(tree.par, phyDat(dat.file), optGamma = T)
            tree.par <- optim.pml(pml.temp, optEdge = T, optInv = T)$tree
        }
        else if (length(grep("+G", best.model.temp)) == 1 && 
            length(grep("+I", best.model.temp)) == 1) {
            pml.temp <- pml(tree.par, phyDat(dat.file), optInv = T, 
                optGamma = T)
            tree.par <- optim.pml(pml.temp, optEdge = T, optGamma = T, 
                optInv = T)$tree
        }
        return(tree.par)
    }
    if (para == T) {
        print("running parallel version, please wait")
        require(doParallel)
        require(foreach)
        print("making clusters")
        cl <- makeCluster(ncore)
        registerDoParallel(cl)
        print("clusters registered")
        i = data.files
        trees.opt <- foreach(dat = data.files, i = 1:length(data.files), 
            .packages = c("phangorn", "ape")) \%dopar\% optim.trees.par(dat)
        stopCluster(cl)
        print("parallelized run complete")
    }
    options(warn = 1)
    for (i in 1:length(trees.opt)) {
        names(trees.opt)[i] <- substr(file.names[i], 1, nchar(file.names[i]) - 
            nchar(form) - 1)
    }
    if (save.trees == T) {
        print("saving trees")
        class(trees.opt) <- "multiPhylo"
		write.tree(trees.opt, file = paste(tree.file.names, ".trees", sep=""), tree.names=T)
    }
    class(trees.opt) <- "multiPhylo"
    l.res <- list()
    l.res[[1]] <- trees.opt
    l.res[[2]] <- model.table
    return(l.res)
  }
}
\keyword{ ML branch length estimation }
\keyword{ phangorn }
