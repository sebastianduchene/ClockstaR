diagnostics.clockstar <- function(groups.obj, save.plots = FALSE, plots.file = "clockstar.diagnostics.pdf", interactive = TRUE){
    if(length(groups.obj) == 3){
        stop("The object specified is not the correct format. Please provide an the output of function get.all.groups. Please use pam.restuls = T")
    }
    gap.stat <- groups.obj[[2]][[1]][, 3]
    ks <- 1:length(gap.stat)
    errs.up <- gap.stat + groups.obj[[2]][[1]][, 4]
    errs.down <- gap.stat - groups.obj[[2]][[1]][, 4]
    k.matrix <- matrix(NA, nrow = nrow(groups.obj[[4]]), ncol = 1)
    for(i in 1:length(ks)){
        part <- pam(groups.obj[[4]], k = ks[i])
        k.matrix <- cbind(k.matrix, as.matrix(part[[3]]))
    }
    k.matrix <- k.matrix[, 2:ncol(k.matrix)]
    plot.matrices <- function(clust.matrix){
        par(mar = c(2, 4, 2, 4))
        image(t(as.matrix(clust.matrix)), axes = F, col = sample(rainbow(ncol(clust.matrix)*100), ncol(clust.matrix)) , ylab = "Data subset")
        mtext(text = rownames(clust.matrix), side = 2, line = 0.3, at = seq(0, 1, 1/(nrow(clust.matrix) - 1)), las = 1, cex = 0.6)
        mtext(text = paste0("k=", ks), side = 3, line = 0.3, at = seq(0, 1, 1/(ncol(clust.matrix) - 1)), las = 1, cex = 0.6)
    }
    if(interactive == TRUE){
    	par(mfrow = c(2, 1))
    	par(mar = c(4, 4, 4, 4))
    	errbar(x = ks, y = gap.stat, yplus = errs.up, yminus = errs.down, type = "b", pch = 20, xlab = expression(italic(k)), ylab = "Gap statistic")
    	points(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
    	lines(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
    	lines(1:(length(groups.obj[[3]]) + 1), y = rep(0, length(groups.obj[[3]]) + 1), col = "blue")
    	legend(x = max(ks)*0.7, y = min(gap.stat)*0.5, legend = c("Gap statistic and S.E." ,"Delta Gap statistic", "Gap statistic = 0"), cex = 0.8, fill = c("black", "red", "blue"), bty = "n")
    	plot.matrices(k.matrix)
    }
    if(save.plots == TRUE){
        pdf(file = plots.file)
        par(mfrow = c(2, 1))
        par(mar = c(4, 4, 4, 4))
        errbar(x = ks, y = gap.stat, yplus = errs.up, yminus = errs.down, type = "b", pch = 20, xlab = expression(italic(k)), ylab = "Gap statistic")
        points(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
        lines(1:length(groups.obj[[3]]), groups.obj[[3]], col = "red", pch = 20)
        lines(1:(length(groups.obj[[3]]) + 1), y = rep(0, length(groups.obj[[3]]) + 1), , col = "blue")
        legend(x = max(ks)*0.7, y = min(gap.stat)*0.5, legend = c("Gap statistic and S.E." ,"Delta Gap statistic", "Gap statistic = 0"), cex = 0.8, fill = c("black", "red", "blue"), bty = "n")
        plot.matrices(k.matrix)
        dev.off()
    }
    return(k.matrix)
}