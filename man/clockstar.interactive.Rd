\name{clockstar.interactive}
\alias{clockstar.interactive}
\title{
Run clockstar interactively
}
\description{
This function can be used to run ClockstaR interactively. It prompts the use for a file with the gene trees in NEWICK format and some options for the clustering methods.
}
\usage{
clockstar.interactive()
}
\author{
Sebastian Duchene
}

\references{
Duchêne, S, Molak, M., and Ho, SYW. "ClockstaR: choosing the number of relaxed-clock models in molecular phylogenetic analysis." Bioinformatics (2013): btt665.
}

\seealso{
optim.trees.interactive() to obtain the gene trees.
}
\examples{
\dontrun{
clockstar.interactive()
}

## The function is currently defined as
function () 
{
    print("Welcome to ClockstaR2")
    trees.file <- readline("please drag or type in the path to your gene trees file in NEWICK format:\n")
    print(paste("reading trees from file: ", trees.file))
    trees <- read.tree(trees.file)
    print(paste("I read", length(trees), "from your file"))
    if (length(names(trees) > 0)) {
        print(paste("The names of your trees are:"))
        print(names(trees))
    }
    else {
        warning("These trees have no names")
    }
    if (("foreach" \%in\% installed.packages()[, 1]) && "doParallel" \%in\% 
        installed.packages()[, 1]) {
        print("Packages foreach and doParallel are available for parallel computation")
        para <- readline("Should we run ClockstaR in parallel (y / n)? (This is good for large data sets)\n")
        if (para == "y") {
            require(foreach)
            require(doParallel)
            ncore <- as.integer(readline("How many cores should ClockstaR use (integer):\n"))
            print(paste("running ClockstaR in parallel with", 
                ncore, "cores"))
            trees.bsd <- bsd.matrix.para(trees, para = T, ncore = ncore)
        }
        else {
            print("Calculating sBSDmin distances between all pairs of trees")
            trees.bsd <- bsd.matrix(trees)
        }
    }
    else {
        print("Calculating sBSDmin distances between all pairs of trees")
        trees.bsd <- bsd.matrix(trees)
    }
    print("I finished calculating the sBSDmin distances between trees\n")
    default.run <- readline("The settings for clustering with ClockstaR are:\nPAM clustering algorithm\nK from 1 to number of data subsets-1\nSEmax criterion to select the optimal k\n500 bootstrap replicates\n Are these correct? (y/n)\n")
    if (default.run == "y") {
        part.data <- partitions(trees.bsd)
    }
    else {
        fun.cluster <- "d"
        while (!(fun.cluster \%in\% c("1", "2", "3"))) {
            fun.cluster <- readline("Please select one of the clustering functions bellow:(1-3)\n(1) PAM\n(2) CLARA\n(3) FANNY\n (See the user manual for package cluster for more details\n")
        }
        if (fun.cluster == "1") {
            fun <- pam
        }
        else if (fun.cluster == "2") {
            fun <- clara
        }
        else if (fun.cluster == "3") {
            fun <- fanny
        }
        max.k <- as.numeric(readline("What should be the maximum k to test (the maximum is the number of data subsets - 1)\n"))
        if (max.k >= ncol(trees.bsd[[2]])) {
            stop("The maximun k should be between 1 and the number of data sets - 1. ABORTING")
        }
        criterion <- readline("Please type in the criterion to select the optimal k (This can be firstSEmax, Tibs2001max, globalSEmax, firstmax, or globalmax)\n")
        while (!(criterion \%in\% c("firstSEmax", "Tibs2001max", 
            "globalSEmax", "firstmax", "globalmax"))) {
            criterion <- readline("You have not selected any of the available criteria. Please select again")
        }
        boot.n <- as.numeric(readline("How many boostrap replicates should be run for the Gap statistic?\n"))
        print("The settings for clustering are complete. The settings are:\n")
        print("cluster function, maximum k, criterion")
        print(c(c("PAM", "CLARA", "FANNY")[as.numeric(fun.cluster)], 
            max.k, criterion))
        part.data <- partitions(trees.bsd, FUN = fun, kmax = max.k, 
            B = boot.n, gap.best = criterion)
    }
    print("ClockstaR has finished running")
    print(paste("The best number of partitions for your data set is:", 
        part.data$best.k))
    save.dat <- readline("Do you wish to save the results in a pdf file?(y/n)\n")
    if (save.dat == "y") {
        fil.name <- readline("What should be the name and path of the output file?\n")
        plot.partitions(part.data, save.plot = T, file.name = fil.name)
    }
    else {
        print("please see the results in the active graphics devices")
        plot.partitions(part.data)
    }
    return(part.data)
  }
}
\keyword{ gene-tree }
\keyword{ molecular-clock }
\keyword{ sBSDmin }
