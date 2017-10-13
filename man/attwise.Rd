\name{attwise}
\alias{attwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function that identifies attribute-wise patterns in a data set
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function identifies unique attribute-wise patterns and their frequencies within each trial as a unit of analysis.
}
\usage{
attwise(df, "participant", "trial", "alternative", "attribute")
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
\item{alternative}{ fixations to different alternatives
%%     ~~Describe \code{x} here~~
}
\item{attribute}{ fixations to diffrent attributes
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

#the same data set (df) as for 'altwise' function should be used

attwisePatterns <- attwise(df, "participant", "trial", "alternative", "attribute")

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~attribute-wise }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~information_search }% __ONLY ONE__ keyword per line
