ClockstaR2
=========

This is a re-write version of ClockstaR. This version should only depend on APE and cluster. Optionally, it can use doParallel and foreach. This does not include clockstar for genomic data sets


What is ClockstaR?
-----------------
ClockstaR is an R package for phylogenetic molecular clock analyses of multi-gene data sets. It uses the patterns of among lineage rate variation for the different genes to select the clock-partitioning strategy. The method uses a phylogenetic tree distance metric and an unsupervised machine learning algorithm to identify the optimal number of clock-partitoins, and which genes should be analysed under each of the partitons. The partitioning stategy selected in ClocsktaR can be used for subsequent molecular clock analysis with programs such as [BEAST](http://beast.bio.ed.ac.uk/Main_Page), [MrBayes](http://mrbayes.sourceforge.net/), [PhyloBayes](http://megasun.bch.umontreal.ca/People/lartillot/www/index.htm),and others.

Please follow [this link](http://bioinformatics.oxfordjournals.org/content/early/2013/12/02/bioinformatics.btt665.full) for the original publication.

This repository contains the program as an R package. It includes tutorials and example data sets. To download follow this [link](https://github.com/sebastianduchene/clockstar) and click on "Download zip" at the bottom right of the page. To install the program follow the steps in Tutorial 1 in the "tutorials" folder. ClockstaR requires an [R](http://www.r-project.org/) installation. It also requires some R dependencies, which can be obtained through R, as explained in Tutorial 1.

Please send any requests or questions to Sebastian Duchene (sebastian.duchene[at]sydney.edu.au). Some other software and resources can be found at the [Molecular Ecology, Evolution, and Phylogenetics Laboratory](http://sydney.edu.au/science/biology/meep/) at the University of Sydney.


Getting started:
----------------

Download a this repository as a zip file and unzip. The following instructions use the example_clockstar_data folder, which contains some fasta files and a tree. These are simulated data under four pattenrns of evolutionary rate variation.

To install directly from github open R and type:

```
install.packages("devtools")
install_github('ClockstaR', 'sebastianduchene')
```

To see an example type:

```
example(ClockstaR2)
```

To optimise trees interactively (drag the clockstar_example_data folder in this repository):

```
optim.trees.interactive()
```

To run clockstar interactively (drag the tree file produced with optim.trees.interactive() ):

```{r}
clockstar.interactive()
```



Pending:
--------

- Update tutorial 
