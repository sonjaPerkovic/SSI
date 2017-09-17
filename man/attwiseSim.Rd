\name{attwiseSim}
\alias{attwiseSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function that a creates random data set and identifies attribute-wise patterns in it
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function creates two new varibles (alterntive and attribute) by sampling evenly from a defined set of numbers and letters. After the function adds these two columns to the original data set, from which the existing alternative and attribute variables have been previously deleted, it identifies attribute-wise patterns and their frequencies within each trial as a unit of analysis.
}
\usage{
attwiseSim(dfRan, "participant", "trial", num_alt, num_att)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfRan}{ object of class data table
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
}
\details{
This function uses the same data set as the function 'attwise'. The difference is that for this function to work, alternative and attribute columns should be deleted from the data set before compiling the function.
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
#the same example as for 'altwiseSim' function can be used
#creating a data set with participant and trial variables only by copying them from df data set
#df data set is created in an example for 'altwise' function
dfRan <- NULL
dfRan$participant <- df$participant
dfRan$trial <- df$trial
dfRan <- as.data.table(dfRan)

attwisePatternsRan <- attwiseSim(dfRan, "participant", "trial", 3, 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~attribute-wise }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~random }% __ONLY ONE__ keyword per line
\keyword{ ~information search }% __ONLY ONE__ keyword per line
