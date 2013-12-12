get.all.groups.k <-
function(data.obj, k.man = 1, save.partitions = F, file.name = "partitions.txt"){
    if(class(data.obj) %in% c("data.frame", "matrix")){
     print("The data object is a data frame or a matrix. The partitions will be defined with the PAM algorithm")
    
    
    if(k.man == 1){
        stop("The selected number of partitions is 1. Please select k.man > 1")
    }else{
    pam.k <- pam(data.obj, k = k.man)
    partitions <- unique(pam.k$clustering)
    res.data <- list()
    for(i in 1:length(partitions)){
        res.data[[i]] <- names(which(pam.k$clustering == partitions[i]))
    }
    names(res.data) <- paste("Partition_",1:length(res.data))
    if(save.partitions==T){
        cat(paste("Partitions selected with manual mode for k =", k.man, "\n"), file=file.name)
        for(m in 1:length(res.data)){
            cat(names(res.data[m]), file=file.name, sep="\n",append=T)
            cat(res.data[[m]], file=file.name,append=T)
            cat("\n", file=file.name,append=T)
        }
    }
    return(res.data)
	}
	}else if(class(data.obj) == "phylo"){
           print("The data object is a dendrogram (phylo object). The partitions will be defined but cutting the dendrobran along the longest edge")
tree <- data.obj
k = 2

tree.list <- list()
tree.list[[1]] <- tree
while(length(tree.list) < k){
	diams <- sapply(tree.list, function(tr){if(class(tr)=="phylo"){ return(max(cophenetic(tr)))}else{return(0)}})
	tree.to.cut <- tree.list[[which(diams==max(diams))]]
	tree.list <- tree.list[-which(diams==max(diams))]
	tree.list[c(length(tree.list)+1,length(tree.list)+2)] <- cut.trees.beta(tree.to.cut, beta=0)
}
for(l in 1:length(tree.list)){
	if(class(tree.list[[l]])=="phylo"){
    	tips <- tree.list[[l]]$tip.label
        tree.list[[l]] <- tips
    }
}
names(tree.list) <- paste("Partition_", 1:length(tree.list))

if(save.partitions==T){
	cat(paste("partitions with selected k =", k, "\n"), file=file.name)
    	for(m in 1:length(tree.list)){
	    	cat(names(tree.list[m]), file=file.name, sep="\n",append=T)
	    	cat(tree.list[[m]], file=file.name,append=T)
	    	cat("\n", file=file.name,append=T)
    	}
    }
    return(tree.list)
}
}
