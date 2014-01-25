min.dist.topo.mat <-
function(tree.list){

	d.mat <- matrix(NA, nrow=length(tree.list), ncol=length(tree.list))
	rownames(d.mat) <- names(tree.list)
	colnames(d.mat) <- names(tree.list)

	s.mat <- d.mat
	print("Estimating tree distances")

	if(length(tree.list) > 3){
		d.mat.lin <- vector()
			d.mat.lin <- sapply(2:nrow(d.mat), function(a){print(paste("estimating distances",a-1, "of", nrow(d.mat)-1)); lapply(tree.list[1:(a-1)], function(y){min.dist.topo(tree1=y, tree2=tree.list[[a]]) })  })

		for(a in 1:length(d.mat.lin)){
			vec.temp.dist <- vector()
			vec.temp.scale <- vector()
			for(b in 1:length(d.mat.lin[[a]])){
				vec.temp.dist[b] <- d.mat.lin[[a]][[b]][1]
				vec.temp.scale[b] <- d.mat.lin[[a]][[b]][2]
			}
			d.mat[a+1,1:length(vec.temp.dist)] <- vec.temp.dist
			s.mat[a+1,1:length(vec.temp.dist)] <- vec.temp.scale
		}


            }else if(length(tree.list) <= 3){
                    stop("The number of gene trees is <= 3. ClockstaR requires at least gene 4 trees")
                }
        d.mat.lin <- min.dist.topo(tree.list[[1]], tree.list[[2]])
        d.mat[2,1] <- d.mat.lin[1]
        s.mat[2,1] <-  d.mat.lin[2]
        res.list <- list()
        res.list[[1]] <- as.dist(d.mat)
        res.list[[2]] <- s.mat
        return(res.list)
}
