\name{altwiseSimRep}
\alias{altwiseSimRep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function that repeats 'altwiseSim' function n times
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function that repeats 'altwiseSim' function n times and row binds created data sets.
}
\usage{
altwiseSimRep(dfRan,"participant", "trial", num_alt, num_att, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfRan}{ random data set created after compiling 'altwiseSim' function
%%     ~~Describe \code{x} here~~
}
\item{participant}{ identifies each unique subject
%%     ~~Describe \code{x} here~~
}
\item{trial}{ identifies each unique trial
%%     ~~Describe \code{x} here~~
}
\item{num_alt}{ number of alternatives
%%     ~~Describe \code{x} here~~
}
\item{num_att}{ number of attributes
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
altwisePatternsSim = altwiseSimRep(dfRan,"participant", "trial", 3, 3, 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~alternative-wise }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~simulation }% __ONLY ONE__ keyword per line
\keyword{ ~information_search }% __ONLY ONE__ keyword per line
