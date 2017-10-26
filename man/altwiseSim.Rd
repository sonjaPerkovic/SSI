\name{altwiseSim}
\alias{altwiseSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function that creates a random data set and identifies alternative-wise patterns in it
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function creates two new varibles (alterntive and attribute) by sampling evenly from a defined set of numbers and letters. After the function replaces   the variables from the original data set with these two newly created alternative and attribute variables, it identifies alternative-wise patterns and their frequencies within each trial as a unit of analysis.
}
\usage{
altwiseSim(df, "participant", "trial", num_alt, num_att)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{ object of class data frame
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
%%  ~~ If necessary, more details than the description above ~~
This function uses the same data set as the function 'altwise'.
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
Sonja Perkovic
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

#data set is the same as in an example for 'altwise' function

set.seed(1)
infoSearch <- data.frame(participant=rep(c(1:50), each = 400), trial=rep(c(1:200), each = 100), alternative = sample(1:4, 20000, T), attribute = sample(c("a","b","c","d"), 20000, T))

altwisePatternsRan = altwiseSim(df, "participant", "trial", 4, 4)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~alternative-wise }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~random }% __ONLY ONE__ keyword per line
\keyword{ ~information_search }% __ONLY ONE__ keyword per line
