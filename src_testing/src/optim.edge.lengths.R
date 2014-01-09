optim.edge.lengths <-
function(directory, fixed.tree ,form="fasta", model.test = F, save.trees=F, tree.file.names="output", para = F, ncore = 1){
	options(warn=-1)
	directory = paste(directory, "/", sep="")
	file.names <- dir(directory)

	drop.tip <- ape::drop.tip

	file.names <- file.names[grepl(form, file.names)]

	model.table <- matrix(NA, nrow=length(file.names), ncol=3)
	colnames(model.table) <- c("file", "BIC", "model")
	model.table[,1] <- file.names
	data.files <- list()
	trees.opt <- list()
	print("reading files")

	for(a in 1:length(file.names)){
		data.files[[a]] <- read.dna(paste(directory, file.names[a], sep=""), format=form)
	}


### Function to populate optim.trees in non-parallel version
if(para == F){
	for(b in 1:length(file.names)){

		tax.keep.temp <- fixed.tree$tip.label %in% rownames(data.files[[b]])
		trees.opt[[b]] <- drop.tip(fixed.tree, as.character(fixed.tree$tip.label[!tax.keep.temp]))
        trees.opt[[b]]$edge.length <- rtree(nrow(data.files[[b]]))$edge.length
		pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]),inv=0, shape=1, k=1)
		print(paste("model testing dataset", file.names[b], b, "of", length(file.names)))
### select to avoid model testing
	if(model.test == T){
		dat.temp <- phyDat(data.files[[b]])
		model.temp <- modelTest(dat.temp, multicore = T)
#		model.temp <-  modelTest(pml.temp, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, multicore = T)


### select to avoid model tetsting
        model.table[b,2:3] <- c(model.temp$BIC[model.temp$BIC==min(model.temp$BIC)], model.temp$Model[model.temp$BIC==min(model.temp$BIC)])
        best.model.temp <- model.temp$Model[model.temp$BIC==min(model.temp$BIC)][1]
		}else if(model.test == F){
			best.model.temp <- "JC" 
		}
		if(length(grep("+G" ,best.model.temp))==0 && length(grep("+I", best.model.temp))==0){
			pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]))
			trees.opt[[b]] <- optim.pml(pml.temp, optEdge=T)$tree
		}else if(length(grep("+G", best.model.temp))==1 && length(grep("+I", best.model.temp))==0){
				pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]), optInv=T)
				trees.opt[[b]] <- optim.pml(pml.temp, optEdge=T, optGamma=T)$tree
			}else if(length(grep("+G", best.model.temp))==0 && length(grep("+I", best.model.temp))==1){
				pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]),optGamma=T )
				trees.opt[[b]] <- optim.pml(pml.temp, optEdge=T, optInv=T)$tree
				}else if(length(grep("+G", best.model.temp))==1 && length(grep("+I", best.model.temp))==1){
					pml.temp <- pml(trees.opt[[b]], phyDat(data.files[[b]]), optInv=T, optGamma=T)
					trees.opt[[b]] <- optim.pml(pml.temp, optEdge=T, optGamma=T, optInv=T)$tree
					}
		print(paste("optimized edge lengths for tree", b ,"of", length(file.names)))

	}
}

###############End non parallel version

################ Begin prallel version 
optim.trees.par <- function(dat){
	tree.par <- fixed.tree
	dat.file <- dat	
	tax.keep.temp <- fixed.tree$tip.label %in% rownames(dat.file)
	tree.par <- drop.tip(fixed.tree, as.character(fixed.tree$tip.label[!tax.keep.temp]))
	lens.temp <- rtree(nrow(dat.file))
	tree.par$edge.length <- lens.temp$edge.length
	pml.temp <- pml(tree.par, phyDat(dat.file), inv=0, shape=1, k=1)

	if(model.test == T){
		dat.temp <- phyDat(dat)
		model.temp <- modelTest(dat.temp, multicore = T)
#		model.temp <-  modelTest(pml.temp, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, multicore = T)
	    best.model.temp <- model.temp$Model[model.temp$BIC==min(model.temp$BIC)][1]
		}else if(model.test == F){
			best.model.temp <- "JC" # disselect to avoid model testing
			}
		if(length(grep("+G" ,best.model.temp))==0 && length(grep("+I", best.model.temp))==0){
			pml.temp <- pml(tree.par, phyDat(dat.file))
			tree.par <- optim.pml(pml.temp, optEdge=T)$tree
		}else if(length(grep("+G", best.model.temp))==1 && length(grep("+I", best.model.temp))==0){
				pml.temp <- pml(tree.par, phyDat(dat.file), optInv=T)
				tree.par <- optim.pml(pml.temp, optEdge=T, optGamma=T)$tree
			}else if(length(grep("+G", best.model.temp))==0 && length(grep("+I", best.model.temp))==1){
				pml.temp <- pml(tree.par, phyDat(dat.file),optGamma=T )
				tree.par <- optim.pml(pml.temp, optEdge=T, optInv=T)$tree
				}else if(length(grep("+G", best.model.temp))==1 && length(grep("+I", best.model.temp))==1){
					pml.temp <- pml(tree.par, phyDat(dat.file), optInv=T, optGamma=T)
					tree.par <- optim.pml(pml.temp, optEdge=T, optGamma=T, optInv=T)$tree
					}
		return(tree.par)
}


################ Start parallel version
if(para == T){
print("running parallel version, please wait")
	require(doParallel)
	require(foreach)
	print("making clusters")
	cl <- makeCluster(ncore)
	registerDoParallel(cl)
	print("clusters registered")
	i = data.files
	trees.opt <- foreach(dat = data.files, i = 1:length(data.files), .packages = c('phangorn', 'ape')) %dopar% optim.trees.par(dat)
	stopCluster(cl)
print("parallelized run complete")
}
############	
	
	
	
	
	options(warn=1)
	
	for(i in 1:length(trees.opt)){
			names(trees.opt)[i] <- substr(file.names[i], 1, nchar(file.names[i])-nchar(form)-1)
	}

	if(save.trees==T){
		print("saving trees")
		class(trees.opt) <- "multiPhylo"
		write.tree(trees.opt, file=paste(tree.file.names,".trees" ,sep=""), tree.names=T)
	}
	class(trees.opt) <- "multiPhylo"
	l.res <- list()
	l.res[[1]] <- trees.opt
	l.res[[2]] <- model.table
	return(l.res)
}
