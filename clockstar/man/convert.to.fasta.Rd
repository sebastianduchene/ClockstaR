\name{convert.to.fasta}
\alias{convert.to.fasta}
\title{
This function is used to convert nexus files to fasta format for clockstar execution
}
\description{
This function is used to convert nexus files to fasta format for clockstar execution. It uses a folder, and reads all the files with names ending in .nex and converts them to fasta format.
}
\usage{
convert.to.fasta(directory)
}
\arguments{
  \item{directory}{
This is the directory where the files should be located
}
}
\details{
This function is used internally by run.clockstar() when the files are in nexus format with .nex termination. In this case specify run.clockstar(format = "nex")
}
\value{
This function does not return a value, instead it converts the nexus files in a directory into fasta format.
}
\references{
Duchene et al. 
}
\author{
Sebastian Duchene
}
\note{
Internal function used by run.clockstar(). Can also be used directly, see example
}

\seealso{
See user manual for ClockstaR
}
\examples{
\dontrun{
data(woodmouse)
system("mkdir clockstar_function_example")
write.nexus.data(as.list(woodmouse), file = "./clockstar_function_example/woodmouse.nex")
convert.to.fasta("./clockstar_function_example")
}
## The function is currently defined as
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
}

\keyword{ fasta format }
\keyword{ nexus format }
