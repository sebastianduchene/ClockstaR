partitions.topomat <-
function(topomat.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", kmax = "",...){
  dimat <- as.matrix(topomat.object[[1]])
  if(kmax == ""){
	kmax <- nrow(dimat) - 1
  }
  dimat.mds <-cmdscale(dimat, k = 3) 
  if(find.best == T){
    clusdat <- clusGap(dimat.mds, B = B, FUNcluster = FUN, K.max = kmax)
    npart <- maxSE(f = clusdat$Tab[, 3], SE.f = clusdat$Tab[, 4], method = gap.best)
  }
  parts.list <- list()
  for(i in 1:kmax){
    clus.temp <- FUN(dimat.mds, k = i, ...)
    parts.list[[i]] <- clus.temp$clustering
  }
  parts.mat <- do.call("cbind", parts.list)
  rownames(parts.mat) <- rownames(dimat.mds)
  colnames(parts.mat) <- paste0("k=", 1:ncol(parts.mat))
  
  res <- list(parts.mat, range(topomat.object[[1]]))
  names(res) <- c("parts.mat", "range.dist") 
  if(exists("npart")){
    colnames(res[[1]])[npart] <- paste0(colnames(parts.mat)[npart], "BEST")
    res[[3]] <- npart 
    names(res)[3] <- "best.k"
  }
  if(find.best == T){
    res[[4]] <- clusdat
    names(res)[4] <- "clusterdata"
  }

  class(res) <- "partitions"
  return(res)
  
}
