run.clockstar <-
function(files.directory = "", out.file.name = "", para = F, ncore = 1, model.test = F, mode.run = "", k.man = 2, k.max,...){
	while(files.directory == ""){	
	#if(files.directory == ""){	
		print("Please select the file with the sequence data and the phylogenetic tree ")
		files.directory <- tk_choose.dir(caption = "Select the folder with the data")
	#}
	}
	if(length(grep("tre|fasta|nex*", dir(files.directory))) == 0){
		stop("There are no .tree, .fasta, or .nexus files in the directory. Please make sure the files have the .tre, .fasta, or .nex ")
	}
		
	if(out.file.name == ""){
		out.file.name <- readline("Please give a name for the results file (default is clockstar.output): ")
		if(out.file.name == ""){
			out.file.name <- "clockstar.output"
		}
	}
	
	while(mode.run == ""){
		mode.run <- readline("Run clockstar in automatic or manual mode (type 'a' or 'm')? ")
	
		if(mode.run == "m"){
			k.man <- readline("Please input a value for k : ")
			k.man <- as.numeric(k.man)
			while(!is.numeric(k.man)){
				k.man <- readline("Please input a NUMERIC value for k : ")	
			}
			print(paste("you have selected a k of: ", k.man))
		}
		if(mode.run == "a"){
			print("ClockstaR will run in automatic mode")
		#	k.max <- readline("Please input a maximum number of partitions (k.max). This should be 1 < k.max < number of data subsets - 1: ")
		#	if(k.max == ""){
		#		k.max <- 2
		#	}
		}
	}
	
 	files <- dir(files.directory)
	files.format <- files[!(grepl(".tre", files) | grepl(out.file.name, files))][1]
	files.format <- substr(files.format, regexpr("[.]", files.format)[1] + 1, nchar(files.format))

	if(files.format == "nex"){
    	convert.to.fasta(files.directory)
	}

        # We now want all the output into a single file. This may only work in unix-like OS
        init.dir <- getwd()
        setwd(files.directory)
        if(any(dir() == out.file.name)){
        	delete.files <- NA
        	while(!(delete.files %in% c("Y", "y", "N", "n"))){
	        	delete.files <- readline("A clockstar run with the same output name already exists in this folder. Overwrite? (y/n) ")
    		}
        	if(delete.files %in% c("Y", "y")){
        		if(Sys.info()[1] == "Windows"){
					shell(paste("rm -rf", out.file.name))
	        		shell(paste("mkdir", out.file.name))        			
        			}else{
					system(paste("rm -rf", out.file.name))
	        		system(paste("mkdir", out.file.name))
	        		}
    		    setwd(out.file.name)        		
        	}else{
        		new.dir <- paste(out.file.name, length(grep(out.file.name, dir())) + 1, sep = "") 
        		if(Sys.info()[1] =="Windows"){
        			shell(paste("mkdir ", new.dir, sep =""))
        			}else{
        				system(paste("mkdir ", new.dir, sep =""))
        			}
        		setwd(new.dir)
        	}
        }else{
        	if(Sys.info()[1] == "Windows"){
        		shell(paste("mkdir", out.file.name))
        		}else{
        			system(paste("mkdir", out.file.name))
        		}
        	setwd(out.file.name)
		}

	fix.tree <- read.tree(paste(files.directory, files[grep(".tre", files)[1]], sep="/"))
	fix.tree$edge.length <- runif(length(fix.tree$tip.label)*2-1)
	print("OPTIMISING BRANCH LENGTHS")
	opt.trees <- optim.edge.lengths(files.directory, fix.tree, save.trees=T, model.test = model.test, tree.file.names="optimized", para = para, ncore = ncore)
	write.table(opt.trees[[2]], file="models.csv" ,sep=",", row.names=F)
	opt.trees.only <- opt.trees[[1]]
	print("FINISHED OPTIMISING BRANCH LENGTHS")
	print("CALCULATING BSD")
	if(para == F){
		bsd.matrix <- min.dist.topo.mat(opt.trees.only)
	}else if(para == T){
		bsd.matrix <- min.dist.topo.mat.para(opt.trees.only, para = T, ncore = ncore)
	}
	print("FINISHED CALCULATING BSD")
	write.table(as.matrix(bsd.matrix[[1]]), file="bsd.distances.csv", sep=",")
	write.table(as.matrix(bsd.matrix[[2]]), file="scaling.factors.csv", sep=",")

	#some data manipulation to bypass errors generated in 2 partition analyses
	if(ncol(bsd.matrix[[2]]) > 2){
		bsd.dendrogram <- nj(bsd.matrix[[1]])
		write.tree(bsd.dendrogram, file="bsd.dendrogram.tre")
	}else if(ncol(bsd.matrix[[2]])==2){
		bsd.matrix[[1]] <- as.matrix(bsd.matrix[[1]])
		bsd.matrix[[1]] <- cbind(c(0.5, 0.5), bsd.matrix[[1]])
		bsd.matrix[[1]] <- rbind(c(0.5,0.5, 0.5), bsd.matrix[[1]])
		colnames(bsd.matrix[[1]])[1] <- 3
		rownames(bsd.matrix[[1]])[1] <- 3
		bsd.matrix[[1]] <- as.dist(bsd.matrix[[1]])
		bsd.dendrogram <- nj(bsd.matrix[[1]])
		bsd.dendrogram <- drop.tip(bsd.dendrogram, "3")
		write.tree(bsd.dendrogram, file="bsd.dendrogram.tre")
	}else{
		print("ERROR IN DENDROGRAM OF TREE DISTANCES")
	}
#######
#######
# mode.run can be "a" or "m"
#######
#######

	bsd.data <- as.matrix(bsd.matrix[[1]])
#	save(groups, file ="grs.rda")
	print("CALCULATING NUMBER OF GROUPS")
	if(mode.run=="a"){
		print(paste("running in automatic mode"))
#		print(getwd())
#		print(mode.run)
#		print(class(bsd.data))

### 		 bsd.matrix is NOT A MATRIX, ITS A LIST

       	print(bsd.data)
       	parts <- get.all.groups(bsd.data, save.partitions=T, pam.results = T, ...)
       	diagnostics.output <- diagnostics.clockstar(parts, save.plots = TRUE, interactive = FALSE, ...)
	}else if(mode.run=="m"){
		print(paste("running manual mode with k =", k.man))
		parts <- get.all.groups.k(bsd.data, k = k.man, save.partitions=T)
	}
	
	
	print("FINISHED CALCULATING NUMBER OF GROUPS")	
	
	print("PLOTTING")
	if(mode.run == "a"){
		parts.pre <- parts
		parts <- parts[[1]]
	}
	#Now plotting all the above

	pdf("bsd.plots.pdf")
	par(mfrow=c(3,1))
# add image plot for the partitions
	hist(bsd.data, xlab=expression(italic(sBSDmin)), main=expression(italic(sBSDmin)), freq=T)

k.expe.vec <- vector()
k.expe.names <- vector()

for(i in 1:length(parts)){
    k.expe.vec <- c(k.expe.vec, parts[[i]])
    k.expe.names <- c(k.expe.names, rep(names(parts)[i],length(parts[[i]])))
}

k.expe.names <- as.numeric(as.factor(k.expe.names))

k.matrix <- k.expe.names
names(k.matrix) <- k.expe.vec

par(mar = c(4, 4, 4, 4))
image(t(as.matrix(k.matrix)), axes = F, col = rainbow(length(k.matrix)), main = "Partition assignments (Colour represents individual partitions)")
mtext(text = names(k.matrix), side = 2, line = 0.3, at = seq(0, 1, 1/(length(k.matrix) - 1)), las = 1, cex = 0.8)
	
	plot(bsd.dendrogram, "unrooted")
	dev.off()
	
	if(mode.run=="m"){
		parts.pre <- get.all.groups(bsd.data, pam.results = TRUE, ...)
		diagnostics.output <- diagnostics.clockstar(parts.pre, save.plots = TRUE, interactive = FALSE, ...)
	}else{
		diagnostics.output <- diagnostics.clockstar(parts.pre, save.plots = TRUE, interactive = FALSE, ...)
	}
	
	setwd(init.dir)
	print("FINISHED RUN")
}
