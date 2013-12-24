get.all.groups <-
function(data.obj, n.b = 500, pam.results = F, save.partitions = F, file.name = "partitions.txt", beta = 0.03, ...){
	if(class(data.obj) %in% c("data.frame", "matrix")){
		kmax <- nrow(data.obj) - 1
    	print("FINDING THE BEST NUMBER OF PARTITIONS (k) WITH GAP STATISTIC AND THE PAM ALGORITHM")
    	print(kmax)
    	pam.gap <- clusGap(data.obj, pam, B = n.b, K.max = kmax)
    	rownames(pam.gap[[1]]) <- paste0("k=", 1:nrow(pam.gap[[1]]))
    	gaps <- pam.gap$Tab[, 3]
    	diffs <- vector()
    	for(i in 1:(length(gaps) - 1)){
        	diffs <- c(diffs, gaps[i] - gaps[i + 1])
    	}
    	names(diffs) <- paste0("k=", 1:length(diffs))
    	all.k <- (1:length(diffs))[diffs > 0]
		best.k <- maxSE(pam.gap[[1]][, 3], pam.gap[[1]][, 4])
    	pam.k <- pam(data.obj, k = best.k, ...)
    	partitions <- unique(pam.k$clustering)
    	res.data <- list()
    	for(i in 1:length(partitions)){
        	res.data[[i]] <- names(which(pam.k$clustering == partitions[i]))
    	}
    	names(res.data) <- paste("Partition_",1:length(res.data))
    	if(save.partitions==T){
        	cat("Partitions selected with automatic mode \n", file=file.name)
        	for(m in 1:length(res.data)){
            	cat(names(res.data[m]), file=file.name, sep="\n",append=T)
            	cat(res.data[[m]], file=file.name,append=T)
            	cat("\n", file=file.name,append=T)
        	}
    	}
    	if(pam.results == T){
        	lis.res <- list(res.data, pam.gap, diffs, data.obj)
        	return(lis.res)
    	}else if(pam.results == F){
        	return(res.data)
        	}
	}else if(class(data.obj) == "phylo"){
    	tree <- data.obj
	    temp.list <- list()
    	min.list <- list()
    	temp.list[[1]] <- tree
    	get.diameter <- function(tr){if(class(tr)=="phylo"){return(max(cophenetic(tr)))}else{return(0) }}
    	diams <- sapply(temp.list, get.diameter)
    	if(any(diams>beta)){
        	while(length(temp.list)!=0){
            	diams <- sapply(temp.list, get.diameter)
            	if(sum(diams <= beta)>0){# If any of the trees have a diameter <= beta then these are saved in min.list
                	diams.beta <- seq(from=1, to=length(diams))[diams <= beta]
                	for(i in diams.beta){
                    	min.list[[length(min.list)+1]] <- temp.list[[i]]
                	}
            	}else{
                	diams.beta=0
            	}
            	if(sum(diams > beta)>0){
                	diams.non.beta <-  seq(from=1, to=length(diams))[diams > beta]
                	temp.list.non.beta <- list()# If any trees have diamters> beta create temp.list.non.beta

                	for(j in 1:length(diams.non.beta)){
                    	temp.list.non.beta[[j]] <- temp.list[[diams.non.beta[j]]]
                	}#save all the trees with diameter > beta to temp.list (rewrite temp.list)
                	temp.list <- temp.list.non.beta

                	sub.list <- list()

                	for(k in 1:length(temp.list)){
                    	cut.temp <- cut.trees.beta(temp.list[[k]], beta)
                    	sub.list[[length(sub.list)+1]] <- cut.temp[[1]]
                    	if(length(cut.temp)==2){
                        	sub.list[[length(sub.list)+1]] <- cut.temp[[2]]
                    	}
                	}
                	temp.list=sub.list[1:length(sub.list)]
            	}else{
                	temp.list <- list()
            	}
        	}
    	}else{
        	min.list <- temp.list
    	}
        	for(l in 1:length(min.list)){
        		if(class(min.list[[l]])=="phylo"){
            		tips <- min.list[[l]]$tip.label
            		min.list[[l]] <- tips
        		}
    		}
    		names(min.list) <- paste("Partition_", 1:length(min.list))

    		if(save.partitions==T){
    			cat(paste("partitions with selected beta =", beta, "\n"), file=file.name)
    			for(m in 1:length(min.list)){
	    			cat(names(min.list[m]), file=file.name, sep="\n",append=T)
	    			cat(min.list[[m]], file=file.name,append=T)
	    			cat("\n", file=file.name,append=T)
    			}
    		}

    	return(min.list)
	}
}
