\name{probAttwise}
\alias{probAttwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function that computes probabilities and probability complements for attribute-wise patterns
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes probabilites and probability complements by comparing the original data set created using 'attwise' function and simulated data sets created using 'attwiseSimRep' function. The prinicple used for comparison is as follows: it is counted how many times does each unique pattern from the original data set occurs that many times within a specific trial in n simulated data sets.
}
\usage{
probAttwise(attwisePatterns, attwisePatternsSim, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{attwisePatterns}{ data set created after applying 'attwise' function
%%     ~~Describe \code{x} here~~
}
\item{attwisePatternsSim}{ data set created after applying 'attwiseSimRep' function n times
%%     ~~Describe \code{x} here~~
}
\item{iter}{ number of iterations
%%     ~~Describe \code{x} here~~
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
Sonja Perkovic
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
function (x)
{
  }
attwiseProbabilities = probAttwise(attwisePatterns, attwisePatternsSim, 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~attribute-wise }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~simulation }% __ONLY ONE__ keyword per line
\keyword{ ~probability }% __ONLY ONE__ keyword per line
\keyword{ ~information_search }% __ONLY ONE__ keyword per line
