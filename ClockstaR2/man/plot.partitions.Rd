\name{plot.partitions}
\alias{plot.partitions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.partitions(partitions.object, save.plot = F, file.name = "results_clockstar")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{partitions.object}{
%%     ~~Describe \code{partitions.object} here~~
}
  \item{save.plot}{
%%     ~~Describe \code{save.plot} here~~
}
  \item{file.name}{
%%     ~~Describe \code{file.name} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (partitions.object, save.plot = F, file.name = "results_clockstar") 
{
    image(t(partitions.object$parts.mat), col = sample(rainbow(ncol(partitions.object$parts.mat))), 
        axes = F, ylab = "Data subset", main = "Partitioning scheeme for values of k")
    mtext(text = rownames(partitions.object$parts.mat), side = 2, 
        line = 0.3, at = seq(0, 1, 1/(nrow(partitions.object$parts.mat) - 
            1)), las = 1, cex = 0.6)
    mtext(text = colnames(partitions.object$parts.mat), side = 3, 
        line = 0.3, at = seq(0, 1, 1/(ncol(partitions.object$parts.mat) - 
            1)), las = 1, cex = 0.6)
    if (save.plot) {
        dev.copy2pdf(file = paste0(file.name, "_matrix.pdf"))
    }
    dev.new()
    plot(partitions.object$clusterdata)
    if (save.plot) {
        dev.copy2pdf(file = paste0(file.name, "_gapstats.pdf"))
        sapply(dev.list(), function(x) dev.off(x))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
