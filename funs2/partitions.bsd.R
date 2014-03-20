

partitions <- function(bsd.object, ...) UseMethod("partitions")



partitions.bsd <- function(bsd.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", ...){
  dimat <- as.matrix(bsd.object[[1]])
  if(!exists("K.max")){
	K.max <- nrow(dimat) - 1
  }

  if(find.best == T){
    clustdat <- clusGap(dimat, B = B, FUNcluster = FUN, K.max = K.max)
    npart <- maxSE(f = clustdat$Tab[, 3], SE.f = clustdat$Tab[, 4], method = gap.best)
  }
  parts.list <- list()
  for(i in 1:(nrow(dimat) - 1)){
    clus.temp <- FUN(dimat, k = i, ...)
    parts.list[[i]] <- clus.temp$clustering
  }
  parts.mat <- do.call("cbind", parts.list)
  rownames(parts.mat) <- rownames(dimat)
  colnames(parts.mat) <- paste0("k=", 1:ncol(parts.mat))
  
  res <- list(parts.mat, range(bsd.object[[1]]))
  names(res) <- c("parts.mat", "range.bsd") 
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
