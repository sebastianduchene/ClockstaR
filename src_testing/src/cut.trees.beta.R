cut.trees.beta <-
function(tree, beta=0.05){
	# Calculate tree diameter
	if(length(tree$tip.label)>=3){
		tree <- unroot(tree)
   }
    tree.diam <- max(cophenetic(tree))
    pruned.trees <- list()

    if(tree.diam>beta & (length(tree$tip.label)>2)){
        di.nodes <- dist.nodes(tree)
        max.edge <- max(tree$edge.length)
        tips <- 1:length(tree$tip.label)
        nodes <- 1:nrow(di.nodes)
        nodes <- nodes[-tips]
        connect.longest.edge <- which(di.nodes==max.edge, arr.ind=T)[1,]

        if((connect.longest.edge[1] %in% nodes) & (connect.longest.edge[2] %in% nodes)){
            taxa.cut1 <- tips(tree, connect.longest.edge[1])
            taxa.cut2 <- tips(tree, connect.longest.edge[2])
            len.tax <- c(length(taxa.cut1), length(taxa.cut2))
            node.cut <- connect.longest.edge[len.tax==min(len.tax)]
            taxa.cut <- tips(tree, node.cut)
        }else if(sum(connect.longest.edge %in% nodes)==1){#If it is a node - tip connection
            node.cut <- connect.longest.edge[connect.longest.edge %in% nodes]
            taxa.cut <- tree$tip.label[connect.longest.edge[connect.longest.edge %in% tips]]
        }

        taxa.prune <- taxa.cut

        # HERE ENDS THE SECTION FOR SELECTING THE NODE AND TAXA TO BE PRUNED
        # After defining the tips to be pruned, create the subtree1. If the number of tips left is >=2
	# then this works just by pruning out the taxa. Otherwise an the object is the tip labels left
        if((length(tree$tip.label)-length(taxa.prune))>=2){
            subtree1 <- drop.tip(tree, taxa.prune)
            if(length(subtree1$tip.label)>=3){
	           	subtree1 <- unroot(subtree1)
  	         }
        }else{
            subtree1 <- tree$tip.label[!(tree$tip.label %in% taxa.prune)]
        }
	# Create the second subtree2 which is the tips that were not included in subtree1
        if(length(taxa.prune) >= 2){
            subtree2 <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% taxa.prune)])
            if(length(subtree2$tip.label)>=3){
            	subtree2 <- unroot(subtree2)
            }
        }else{
            subtree2 <- taxa.prune
        }

        pruned.trees[[1]] <- subtree1
        pruned.trees[[2]] <- subtree2

    }else if(tree.diam < beta){
        pruned.trees[[1]] <- tree
    }else{
    	pruned.trees[[1]] <- tree$tip.label[1]
    	pruned.trees[[2]] <- tree$tip.label[2]
    }
    return(pruned.trees)
}
