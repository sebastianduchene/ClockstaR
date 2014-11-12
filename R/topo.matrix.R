topo.matrix <-
function(tree.list){

  topo_dist_mat <- matrix(NA, length(tree.list), length(tree.list))
  rownames(topo_dist_mat) <- gsub('[.].+$', '', names(tree.list))
  colnames(topo_dist_mat) <- rownames(topo_dist_mat)

  print('Estimating tree topology distances')

  if(length(tree.list) > 3) {

    for(i in 2:nrow(topo_dist_mat)){
      print(paste('Estimating topology distances for tree', rownames(topo_dist_mat)[i], i-1, 'of', nrow(topo_dist_mat) - 1))
      temp_trees_dist <- sapply(tree.list[1:i-1], function(x) dist.topo(tree.list[[i]], x))
      topo_dist_mat[i, 1:length(temp_trees_dist)] <- temp_trees_dist
    }

  }else if(length(tree.list) <= 3){
      stop('The number of gene trees is <=3. ClockstaR requires at least 4 gene trees')
  }

  res.data <- list(topo.mat = as.dist(topo_dist_mat))

  class(res.data) <- 'topomat'

  return(res.data)
}
