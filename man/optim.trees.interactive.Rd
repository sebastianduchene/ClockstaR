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
The function does require any arguments. When it is called, it prompts the user for the folder with the data and the substitution models for each data set. In some versions of R it may be necessary to specify the folder with the fasta files and the tree topology. This can be done by using the arguments folder.parts. For example optim.trees.interactive(folder.parts = 'clockstar_example_data')
}
\author{
Sebastian Duchene
}

\examples{
\dontrun{
optim.trees.interactive() # Follow the instructions in the prompt
}

}

