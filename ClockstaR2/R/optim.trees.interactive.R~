optim.trees.interactive <-
function(){

if(!("phangorn" %in% installed.packages()[,1])){
  stop("Package phangorn is not installed. Please install it to run this function")
}

require(phangorn)
folder.parts <- readline("Please drag a folder with the data subsets and a tree topology. The files should be in FASTA format, and the trees in NEWICK\n")
dir.init <- getwd()
setwd(folder.parts)

files.parts <- grep("[.]fas", dir(), value = T)

tree.parts <- grep("[.]tre", dir(), value = T)

if(length(files.parts) < 1 || length(tree.parts) < 1){
  stop("There are no fasta or tree files in this folder. Please make sure that they have the extensions .fasta or .tree")
}

print(paste("I have read", length(files.parts), "data files:\n"))
print(files.parts)

print(paste("I am using the tree in file", tree.parts[1], "to optimise the branch lengths"))

data.files <- list()

for(i in 1:length(files.parts)){
  data.files[[i]] <- read.phyDat(files.parts[i], format = "fasta")
}

fix.tree <- read.tree(tree.parts[1])

out.file.name <- readline("What should be the nane of the file to save the optimised trees?\n")

print("I have read all the files")

print("Please select the substitution model for each data subset")

opt.list <- list()

for(k in 1:length(data.files)){
      pml.temp <- pml(fix.tree, data.files[[i]])
      print("The available substitution models are: JC, F81, K80, HKY, TrNe, TrN, TPM1, K81, TPM1u, TPM2, TPM2u, TPM3, TPM3u, TIM1e, TIM1, TIM2e, TIM2, TIM3e, TIM3, TVMe, TVM, SYM and GTR")
      mod <- readline(paste("Please type in one of the available substitution models", "for data set", files.parts[k], "\n"))
      o.gam <- as.logical(readline("Do you wish to optimise gamma? (T/F)\n"))
      o.k <- as.numeric(readline("How many intervals should be used for the gamma distribution (only applicable when gamma is T) (integer)\n")) 
      o.inv <- as.logical(readline("Do you wish to optimise I (T/F)\n"))
      if(any(c(mod, o.gam, o.k, o.inv) == "")){
        setwd(dir.init)
        stop("Some of the model parameters are empty. ABORTING")
      }
      opt.list[[k]] <- optim.pml(pml.temp, optInv = o.inv, k = o.k, optGamma = o.gam, model = mod)$tree
      plot(opt.list[[k]])
}   

class(opt.list) <- "multiPhylo"

names(opt.list) <- files.parts

write.tree(opt.list, file = out.file.name, tree.names = T)

print(paste("The trees have been saved in ", out.file.name, "in the working directory"))

setwd(dir.init)


}
