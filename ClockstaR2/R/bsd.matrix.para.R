bsd.matrix.para <-
function(tree.list, para = F, ncore = 1){
  require(foreach)
  require(doParallel)
  
    if(length(tree.list) <= 3){
        stop("The number of gene trees is <= 3. ClockstaR requires at least gene 4 trees")
    }
######################
######################
######################
######################
bsd.dist <- function(tree1 , tree2){

	list.tr <- list()
	list.tr[[1]] <- tree1
	list.tr[[2]] <- tree2
        # In every case the tree that is rescaled is tree2, which is the shortest
	lens <- c(sum(tree1$edge.length), sum(tree2$edge.length))
	tree1 <- list.tr[lens==max(lens)][[1]]
	tree2 <- list.tr[lens==min(lens)][[1]]

        # This is the objective function for the optimisation implementation with optim()
        tree.dist.opt <- function(x){
            tree3 <- tree2
            tree3$edge.length <- tree2$edge.length*x
            return(dist.topo(tree1, tree3, method="score"))
        }

        # Optimisation of the tree distance with a starting value of 0
        opt.dist <- optim(0, fn=tree.dist.opt, method = "Brent", lower = 0, upper = 50)

        min.bdi <- opt.dist$value
        scaling <- opt.dist$par

        # Scaling for tree2
        tree2.scaled <- tree2
        tree2.scaled$edge.length <- tree2$edge.length * scaling

        # The trees are scaled so that the mean branch length is 0.05
        # This is an arbitrary value, but is useful so that the total tree lengths are
        # in similar scales
	root.scaling <- 0.05 / mean(c(tree1$edge.length[tree1$edge.length > 0.00001] , tree2.scaled$edge.length[tree2.scaled$edge.length > 0.00001]))

	tree1.root.scaled <- tree1
	tree2.root.scaled <- tree2.scaled

	tree1.root.scaled$edge.length <- tree1$edge.length * root.scaling
	tree2.root.scaled$edge.length <- tree2.scaled$edge.length * root.scaling

	min.bdi.root.scaled <- dist.topo(tree1.root.scaled, tree2.root.scaled, method="score")
        res.vect <- c(min.bdi.root.scaled, scaling, min.bdi)
        names(res.vect) <- c("min.bdi.scaled", "scaling.factor", "min.bdi")
        return(res.vect)
}

	sub.trees <- list()
	for(k in 2:length(tree.list)){
		sub.trees[[k]] <- tree.list[1:k-1]
	}

		compute.tree.dists <- function(tree.sub.list, fix.tree){
			res <- sapply(tree.sub.list, function(a){ return(bsd.dist(fix.tree, a))})
			return(res)
		}

	if(para == T){
		cl <- makeCluster(ncore)
		registerDoParallel(cl)
		print(paste("Clusters registered as follows: ", cl))
		res.par <- foreach(s.trees = sub.trees, j = 1:length(tree.list), .packages = "ape") %dopar% compute.tree.dists(tree.sub.list = s.trees, fix.tree = tree.list[[j]])
		stopCluster(cl)
		}else if(para == F){
			res.par <- foreach(s.trees = sub.trees, j = 1:length(tree.list), .packages = "ape") %do% compute.tree.dists(tree.sub.list = s.trees, tree.list[[j]])
		}



	res.list <- list()

	res.list[[1]] <- matrix(NA, nrow = length(tree.list)  , ncol = length(tree.list) )
	for(m in 2:nrow(res.list[[1]])){
		res.list[[1]][m, 1:ncol(res.par[[m]])] <- res.par[[m]][1,]
	}
	rownames(res.list[[1]]) <- names(tree.list)
	colnames(res.list[[1]]) <- names(tree.list)
	res.list[[1]] <- as.dist(res.list[[1]])

	res.list[[2]] <- matrix(NA, nrow = length(tree.list)  , ncol = length(tree.list) )
	for(m in 2:nrow(res.list[[2]])){
		res.list[[2]][m, 1:ncol(res.par[[m]])] <- res.par[[m]][2,]
	}
	rownames(res.list[[2]]) <- names(tree.list)
	colnames(res.list[[2]]) <- names(tree.list)
	class(res.list) <- "bsd"
	return(res.list)

}
