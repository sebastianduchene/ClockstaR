\name{run.clockstar}
\alias{run.clockstar}
\title{
This function is wrapper for the individual ClockstaR functions. It uses a folder that contains the individual gene alignments and the phylogenetic tree to calculate the BSDmin distance. }
\description{
The folder is assigned through a pop up window. It should contain the individual gene alignments in fasta or nexus format and a phylogenetic tree for the relationships among the taxa. The taxon names should be the same in all alignments and in the tree.
The output is a folder with a pdf file with some useful plots, two matrices for the pairwise BSDmin distances and the scaling factors for the trees, a dendrogram for the BSDmin distances among genes, the substitution models selected, and the individual gene trees.
}
\usage{
run.clockstar(files.directory = "", out.file.name = "", para = F, ncore = 1, model.test = F, mode.run = "", k.man = 2, k.max, ...)
}
\arguments{
  \item{files.directory}{
	This is the name of the input folder. This folder should contain the sequence alignment files in fasta or nexus format, and a tree file in newick format.
}
  \item{out.file.name}{
	This is the name of the output file, where the results will be stored.
}
  \item{para}{
	Use para = T to run the program in parallel. Packages foreach and doParallel are required. When run in parallel with model.test = T clockstar will internally test the substitution model and estimate the branch lengths accordingly. However, a table for the models and BIC scores is not printed.
}
  \item{ncore}{
	ncore is used to select the number of cores when the program is run in parallel.
}
  \item{model.test}{
	To test the substitution model for each alignment select model.test = T. Default is model.test = F.
}
  \item{mode.run}{
	This option allows the user to select the number of groups 'a priori' using mode.run = "k", or using a threshold for the BSDmin distance among trees with mode.run = "beta".
}
  \item{k.man}{
	If mode.run == "m" select the number of groups. The default is k.man = 2
}
  \item{k.max}{
	If mode.run == "a" select the maximum number of groups. The default is the number of data subsets - 1.
}
  \item{\dots}{
	Other arguments passed to other functions
}
}
\details{
This function can be called with no arguments as run.clockstar(). In this case the user is prompted to select the input and output folders, and the values for k and beta. 
}
\value{
This function outputs a folder with several files:
bsd.plots.pdf has plots for the BSDmin histogram, a coloured plot representing the partitions, the sBSDmin dendrogram.
}
\references{
Duchene et al. ClockstaR: Choosing the number of relaxed-clock models in molecular phylogenetic analysis
}
\author{
Sebastian Duchene
}
\note{
contact at clockstar website
}

\seealso{
See user manual available at: clockstar website
}
\examples{
###############
\dontrun{
# create folder to save some simulated data and set the working directory
initial.wd <- getwd()
system("mkdir clockstar_example")
setwd("clockstar_example")

#simulate trees with three different patterns of among-lineage rate variation
# see the example for the function get.all.groups for more details on these simulations.
fixed.topology <- rtree(10)
fixed.topology$edge.length <- NULL
tree.list <- list()
set.seed(1234)
patterns1 <- abs(rnorm(18, 0, 0.5))
set.seed(3456)
patterns2 <- abs(rnorm(18, 0, 0.5))
set.seed(0987)
patterns3 <- abs(rnorm(18, 0, 0.5))
for(i in 1:3){
	tree.list[[i]] <- fixed.topology
	tree.list[[i]]$edge.length <- abs(patterns1 + rnorm(18, 0, 0.05))
}

for(i in 4:7){
	tree.list[[i]] <- fixed.topology
	tree.list[[i]]$edge.length <- abs(patterns2 + rnorm(18, 0, 0.05))
}

for(i in 8:10){
	tree.list[[i]] <- fixed.topology
	tree.list[[i]]$edge.length <- abs(patterns3 + rnorm(18, 0, 0.05))
}


# Simulate sequence data along each gene tree
secs.data <- list()
for(i in 1:10){
    secs.data[[i]] <- as.DNAbin(simSeq(tree.list[[i]], l = 2000))
}

# Write tree and data to the folder created above
write.tree(fixed.topology, file = "simulation.tre")

for(i in 1:10){
    write.dna(secs.data[[i]], format = "fasta", nbcol = -1, colsep = "", file = paste("secs", i, ".fasta", sep = ""))
}
##############
# Run clockstar using the simulated data in manual mode
#run.clockstar.mod(files.directory = "", out.file.name = "", para = F, ncore = 1, model.test = F, mode.run = "", k.man = 2, k.max, ...)

run.clockstar(files.directory = getwd(), out.file.name = "clockstar_test_out", mode.run = "a")

# See the results in the "clockstar_test_out" folder
setwd(initial.wd)

}

## The function is currently defined as
run.clockstar <-
function(files.directory = "", out.file.name = "", para = F, ncore = 1, model.test = F, mode.run = "", k.man = 2, k.max, ...){
	while(files.directory == ""){	
		print("Please select the file with the sequence data and the phylogenetic tree ")
		files.directory <- tk_choose.dir(caption = "Select the folder with the data")
	}
	if	(length(grep("tre|fasta|nex*", dir(files.directory))) == 0){
		stop("There are no tree, fasta, or nexus files in the directory. Please make sure the files have the .tre, .fasta, or .nex ")
	}
	if	(out.file.name == ""){
		out.file.name <- readline("Please give a name for the results file (default is clockstar.output): ")
		if(out.file.name == ""){
			out.file.name <- "clockstar.output"
		}
	}
    while (mode.run == "") {
        mode.run <- readline("Run clockstar in automatic or manual mode (type 'a' or 'm')? ")
        if (mode.run == "m") {
            k.man <- readline("Please input a value for k : ")
            k.man <- as.numeric(k.man)
            while (!is.numeric(k.man)) {
                k.man <- readline("Please input a NUMERIC value for k : ")
            }
            print(paste("you have selected a k of: ", k.man))
        }
        if (mode.run == "a") {
            print("ClockstaR will run in automatic mode")
        }
    }
    files <- dir(files.directory)
    files.format <- files[!(grepl(".tre", files) | grepl(out.file.name, 
        files))][1]
    files.format <- substr(files.format, regexpr("[.]", files.format)[1] + 
        1, nchar(files.format))
    if (files.format == "nex") {
        convert.to.fasta(files.directory)
    }
    init.dir <- getwd()
    setwd(files.directory)
    if (any(dir() == out.file.name)) {
        delete.files <- NA
        while (!(delete.files \%in\% c("Y", "y", "N", "n"))) {
            delete.files <- readline("A clockstar run with the same output name already exists in this folder. Overwrite? (y/n) ")
        }
        if (delete.files \%in\% c("Y", "y")) {
            if (Sys.info()[1] == "Windows") {
                shell(paste("rm -rf", out.file.name))
                shell(paste("mkdir", out.file.name))
            }
            else {
                system(paste("rm -rf", out.file.name))
                system(paste("mkdir", out.file.name))
            }
            setwd(out.file.name)
        }
        else {
            new.dir <- paste(out.file.name, length(grep(out.file.name, 
                dir())) + 1, sep = "")
            if (Sys.info()[1] == "Windows") {
                shell(paste("mkdir ", new.dir, sep = ""))
            }
            else {
                system(paste("mkdir ", new.dir, sep = ""))
            }
            setwd(new.dir)
        }
    }
    else {
        if (Sys.info()[1] == "Windows") {
            shell(paste("mkdir", out.file.name))
        }
        else {
            system(paste("mkdir", out.file.name))
        }
        setwd(out.file.name)
    }
    fix.tree <- read.tree(paste(files.directory, files[grep(".tre", 
        files)[1]], sep = "/"))
    fix.tree$edge.length <- runif(length(fix.tree$tip.label) * 
        2 - 1)
    print("OPTIMISING BRANCH LENGTHS")
    opt.trees <- optim.edge.lengths(files.directory, fix.tree, 
        save.trees = TRUE, model.test = model.test, tree.file.names = "optimized", 
        para = para, ncore = ncore)
    write.table(opt.trees[[2]], file = "models.csv", sep = ",", 
        row.names = FALSE)
    opt.trees.only <- opt.trees[[1]]
    print("FINISHED OPTIMISING BRANCH LENGTHS")
    print("CALCULATING BSD")
    if (para == F) {
        bsd.matrix <- min.dist.topo.mat(opt.trees.only)
    }
    else if (para == T) {
        bsd.matrix <- min.dist.topo.mat.para(opt.trees.only, 
            para = T, ncore = ncore)
    }
    print("FINISHED CALCULATING BSD")
    write.table(as.matrix(bsd.matrix[[1]]), file = "bsd.distances.csv", 
        sep = ",")
    write.table(as.matrix(bsd.matrix[[2]]), file = "scaling.factors.csv", 
        sep = ",")
    if (ncol(bsd.matrix[[2]]) > 2) {
        bsd.dendrogram <- nj(bsd.matrix[[1]])
        write.tree(bsd.dendrogram, file = "bsd.dendrogram.tre")
    }
    else if (ncol(bsd.matrix[[2]]) == 2) {
        bsd.matrix[[1]] <- as.matrix(bsd.matrix[[1]])
        bsd.matrix[[1]] <- cbind(c(0.5, 0.5), bsd.matrix[[1]])
        bsd.matrix[[1]] <- rbind(c(0.5, 0.5, 0.5), bsd.matrix[[1]])
        colnames(bsd.matrix[[1]])[1] <- 3
        rownames(bsd.matrix[[1]])[1] <- 3
        bsd.matrix[[1]] <- as.dist(bsd.matrix[[1]])
        bsd.dendrogram <- nj(bsd.matrix[[1]])
        bsd.dendrogram <- drop.tip(bsd.dendrogram, "3")
        write.tree(bsd.dendrogram, file = "bsd.dendrogram.tre")
    }
    else {
        print("ERROR IN DENDROGRAM OF TREE DISTANCES")
    }
    bsd.data <- as.matrix(bsd.matrix[[1]])
    print("CALCULATING NUMBER OF GROUPS")
    if (mode.run == "a") {
        print(paste("running in automatic mode"))
        print(bsd.data)
        parts <- get.all.groups(bsd.data, save.partitions = T, ...)
        diagnostics.output <- diagnostics.clockstar(parts, save.plots = TRUE, interactive = FALSE, ...)
    }
    else if (mode.run == "m") {
        print(paste("running manual mode with k =", k.man))
        parts <- get.all.groups.k(bsd.data, k = k.man, save.partitions = T)
    }
    print("FINISHED CALCULATING NUMBER OF GROUPS")
    print("PLOTTING")
    if (mode.run == "a"){
		parts.pre <- parts
		parts <- parts[[1]]
	}
    pdf("bsd.plots.pdf")
    par(mfrow = c(3, 1))
    hist(bsd.data, xlab = expression(italic(sBSDmin)), main = expression(italic(sBSDmin)), 
        freq = T)
    k.expe.vec <- vector()
    k.expe.names <- vector()
    for (i in 1:length(parts)) {
        k.expe.vec <- c(k.expe.vec, parts[[i]])
        k.expe.names <- c(k.expe.names, rep(names(parts)[i], 
            length(parts[[i]])))
    }
    k.expe.names <- as.numeric(as.factor(k.expe.names))
    k.matrix <- k.expe.names
    names(k.matrix) <- k.expe.vec
    image(t(as.matrix(k.matrix)), axes = F, col = rainbow(length(k.matrix)), 
        main = "Partition assignments (Colour represents individual partitions)")
    mtext(text = names(k.matrix), side = 2, line = 0.3, at = seq(0, 
        1, 1/(length(k.matrix) - 1)), las = 1, cex = 0.8)
    plot(bsd.dendrogram, "unrooted")
    dev.off()
    
    if (mode.run=="m") {
		parts.pre <- get.all.groups(bsd.data, pam.results = TRUE, ...)
		diagnostics.output <- diagnostics.clockstar(parts.pre, save.plots = TRUE, interactive = FALSE, ...)
	} else {
		diagnostics.output <- diagnostics.clockstar(parts.pre, save.plots = TRUE, interactive = FALSE, ...)
	}
    
    setwd(init.dir)
    print("FINISHED RUN")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
