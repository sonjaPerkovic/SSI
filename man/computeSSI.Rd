\name{computeSSI}
\alias{computeSSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function that computes SSI
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function takes information from the data sets created with 'probAltwise' and 'probAttwise' functions and computes SSI. It also uses the original data file to compute the string length for each trial.
}
\usage{
computeSSI(df, altwiseProbabilities, attwiseProbabilities, "participant", "trial", "alternative", "attribute")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{ object of class data table
%%     ~~Describe \code{x} here~~
}
\item{altwiseProbabilities}{ data set created using 'probAltwise' function
%%     ~~Describe \code{x} here~~
}
\item{attwiseProbabilities}{ data set created using 'probAttwise' function
%%     ~~Describe \code{x} here~~
}
\item{participant}{ identifies each unique subject
%%     ~~Describe \code{x} here~~
}
\item{trial}{ identifies each unique trial
%%     ~~Describe \code{x} here~~
}
\item{alternative}{ fixations to different alternatives
%%     ~~Describe \code{x} here~~
}
\item{attribute}{ fixations to different attributes
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
#for this function to work, all the previously presented functions need to be compiled
SSI_output = SSI(df, altwiseProbabilities, attwiseProbabilities, "participant", "trial", "alternative", "attribute")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~alternative-wise }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~attribute-wise }% __ONLY ONE__ keyword per line
\keyword{ ~SSI }% __ONLY ONE__ keyword per line
\keyword{ ~information_search }% __ONLY ONE__ keyword per line
