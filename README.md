ClockstaR
=========
Last update: 18 December 2013

What is ClockstaR?
-----------------
ClockstaR is an R package for phylogenetic molecular clock analyses of multi-gene data sets. It uses the patterns of among lineage rate variation for the different genes to select the clock-partitioning strategy. The method uses a phylogenetic tree distance metric and an unsupervised machine learning algorithm to identify the optimal number of clock-partitoins, and which genes should be analysed under each of the partitons. The partitioning stategy selected in ClocsktaR can be used for subsequent molecular clock analysis with programs such as [BEAST](http://beast.bio.ed.ac.uk/Main_Page), [MrBayes](http://mrbayes.sourceforge.net/), [PhyloBayes](http://megasun.bch.umontreal.ca/People/lartillot/www/index.htm),and others.

Please follow [this link](http://bioinformatics.oxfordjournals.org/content/early/2013/12/02/bioinformatics.btt665.full) for the original publication.

This repository contains the program as an R package. It includes tutorials and example data sets. To download follow this [link](https://github.com/sebastianduchene/clockstar) and click on "Download zip" at the bottom right of the page. To install the program follow the steps in Tutorial 1 in the "tutorials" folder. ClockstaR requires an [R](http://www.r-project.org/) installation. It also requires some R dependencies, which can be obtained through R, as explained in Tutorial 1.

Please send any requests or questions to Sebastian Duchene (sebastian.duchene[at]sydney.edu.au). Some other software and resources can be found at the website for the [Molecular Ecology, Evolution, and Phylogenetics Laboratory](http://sydney.edu.au/science/biology/meep/) at the University of Sydney.

News and latest versions
------------------------


Requested features for ClockstaR
--------------------------------
These are features that users have requested. We will make an effort to include them in future versions of the program.

- 18 December 2013
  - Allow the user to select substitution models.
  - Include a tutorial so that the user can input a list of tree with branch lengths estimated with an other program, such as [RAxML](http://www.exelixis-lab.org/)
  - The program should print an error when there are only two data subsets in the analysis. Instead, the current version always returns that the data should be analysed under a single clock-partition. This due to the implementation of the Gap statistic and PAM algorithm used to find the partitoins.
