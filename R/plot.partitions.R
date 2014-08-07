plot.partitions <-
function(partitions.object, save.plot = F, file.name = "results_clockstar"){
  image(t(partitions.object$parts.mat), col = sample(rainbow(ncol(partitions.object$parts.mat))), axes = F, ylab = "Data subset", main = "Partitioning scheme for values of k")
  mtext(text = gsub('([.]([a-z]|[A-Z])+)$', '', rownames(partitions.object$parts.mat)), side = 2, line = 0.3, at = seq(0, 1, 1/(nrow(partitions.object$parts.mat) - 1)), las = 1, cex = 0.6)
  mtext(text = colnames(partitions.object$parts.mat), side = 3, line = 0.3, at = seq(0, 1, 1/(ncol(partitions.object$parts.mat) - 1)), las = 1, cex = 0.6)
  if(save.plot){
	dev.copy2pdf(file = paste0(file.name, "_matrix.pdf"))
  }
  dev.new()
  plot(partitions.object$clusterdata)
  if(save.plot){
	dev.copy2pdf(file = paste0(file.name, "_gapstats.pdf"))
	sapply(dev.list(), function(x) dev.off(x))
  }

}
