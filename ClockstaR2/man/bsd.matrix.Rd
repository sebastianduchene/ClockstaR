\name{bsd.matrix}
\alias{bsd.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bsd.matrix(tree.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree.list}{
%%     ~~Describe \code{tree.list} here~~
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
function (tree.list) 
{
    d.mat <- matrix(NA, nrow = length(tree.list), ncol = length(tree.list))
    rownames(d.mat) <- names(tree.list)
    colnames(d.mat) <- names(tree.list)
    s.mat <- d.mat
    print("Estimating tree distances")
    if (length(tree.list) > 3) {
        d.mat.lin <- vector()
        d.mat.lin <- sapply(2:nrow(d.mat), function(a) {
            print(paste("estimating distances", a - 1, "of", 
                nrow(d.mat) - 1))
            lapply(tree.list[1:(a - 1)], function(y) {
                bsd.dist(tree1 = y, tree2 = tree.list[[a]])
            })
        })
        for (a in 1:length(d.mat.lin)) {
            vec.temp.dist <- vector()
            vec.temp.scale <- vector()
            for (b in 1:length(d.mat.lin[[a]])) {
                vec.temp.dist[b] <- d.mat.lin[[a]][[b]][1]
                vec.temp.scale[b] <- d.mat.lin[[a]][[b]][2]
            }
            d.mat[a + 1, 1:length(vec.temp.dist)] <- vec.temp.dist
            s.mat[a + 1, 1:length(vec.temp.dist)] <- vec.temp.scale
        }
    }
    else if (length(tree.list) <= 3) {
        stop("The number of gene trees is <= 3. ClockstaR requires at least gene 4 trees")
    }
    d.mat.lin <- bsd.dist(tree.list[[1]], tree.list[[2]])
    d.mat[2, 1] <- d.mat.lin[1]
    s.mat[2, 1] <- d.mat.lin[2]
    res.list <- list()
    res.list[[1]] <- as.dist(d.mat)
    res.list[[2]] <- s.mat
    class(res.list) <- "bsd"
    return(res.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
