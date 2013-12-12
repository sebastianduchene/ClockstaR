convert.to.fasta <-
function(directory){
	d <- getwd()
	setwd(directory)
	files <- grep(".nex", dir( ), value=T)
	for(i in 1:length(files)){
	file.temp <- read.nexus.data(files[i])
	write.dna(as.DNAbin(file.temp),file=paste(substr(files[i], 1,nchar(files[i])-4), ".fasta", sep=""),  format="fasta", nbcol=-1, colsep="")
	system(paste("rm", files[i]))
	}
	setwd(d)
}
