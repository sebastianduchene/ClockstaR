\name{optim.trees.interactive}
\alias{optim.trees.interactive}
\title{
Optimise the branch lengths of gene trees interactively
}
\description{
optim.tree.interactive can be used to optimise the gene trees interactively. The data should be a folder with the alignment for each data subset (or gene) and a tree topology. The alignments should be in FASTA format and the tree should be in NEWICK. Please see the example data set.
}
\usage{
optim.trees.interactive()
}
\details{
The function does not take any arguments. When it is called, it prompts the user for the folder with the data and the substitution models for each data set.
}
\author{
Sebastian Duchene
}

\examples{
\dontrun{
optim.trees.interactive() # Follow the instructions in the prompt
}

}

