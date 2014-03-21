\name{partitions.bsd}
\alias{partitions.bsd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
partitions.bsd(bsd.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", kmax = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bsd.object}{
%%     ~~Describe \code{bsd.object} here~~
}
  \item{FUN}{
%%     ~~Describe \code{FUN} here~~
}
  \item{find.best}{
%%     ~~Describe \code{find.best} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{gap.best}{
%%     ~~Describe \code{gap.best} here~~
}
  \item{kmax}{
%%     ~~Describe \code{kmax} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (bsd.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", 
    kmax = "", ...) 
{
    dimat <- as.matrix(bsd.object[[1]])
    if (kmax == "") {
        kmax <- nrow(dimat) - 1
    }
    if (find.best == T) {
        clusdat <- clusGap(dimat, B = B, FUNcluster = FUN, K.max = kmax)
        npart <- maxSE(f = clusdat$Tab[, 3], SE.f = clusdat$Tab[, 
            4], method = gap.best)
    }
    parts.list <- list()
    for (i in 1:kmax) {
        clus.temp <- FUN(dimat, k = i, ...)
        parts.list[[i]] <- clus.temp$clustering
    }
    parts.mat <- do.call("cbind", parts.list)
    rownames(parts.mat) <- rownames(dimat)
    colnames(parts.mat) <- paste0("k=", 1:ncol(parts.mat))
    res <- list(parts.mat, range(bsd.object[[1]]))
    names(res) <- c("parts.mat", "range.bsd")
    if (exists("npart")) {
        colnames(res[[1]])[npart] <- paste0(colnames(parts.mat)[npart], 
            "BEST")
        res[[3]] <- npart
        names(res)[3] <- "best.k"
    }
    if (find.best == T) {
        res[[4]] <- clusdat
        names(res)[4] <- "clusterdata"
    }
    class(res) <- "partitions"
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
