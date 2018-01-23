\name{computeSSI}
\alias{computeSSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function that computes SSI
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function combines integrates all the necessary functions and computes the SSI. It also creates a file with the SSI values for all the alternative- and attribute-wise patterns found. It is important that the format of a data set corresponds to the one provided in the example. Also, column names should match the column names provided in the example. One can adjust the number of alternatives and attributes, number of simulations as well as the minimum threshold for pattern length by setting the function parameters.
}
\usage{
computeSSI(df, dfRan, participant, trial, alternative, attribute, num_alt, num_att, threshold, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{ object of class data frame
%%     ~~Describe \code{x} here~~
}
\item{dfRan}{ the same object which is used for creating random data
%%     ~~Describe \code{x} here~~
}
\item{participant}{ identifies each unique subject
%%     ~~Describe \code{x} here~~
}
\item{trial}{ identifies each unique trial
%%     ~~Describe \code{x} here~~
}
\item{alternative}{ column with eye fixations to different alternatives
%%     ~~Describe \code{x} here~~
}
\item{attribute}{ column with eye fixations to different attributes
%%     ~~Describe \code{x} here~~
}
\item{num_alt}{ number of alternatives in the experiment
%%     ~~Describe \code{x} here~~
}
\item{num_att}{ number of attributes in the experiment
%%     ~~Describe \code{x} here~~
}
\item{iter}{ number of simulation iterations
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
#IMPORTANT! variables representing participants, trials, alternatives and attributes in your data should have the names that match the ones in the example below!

df = data.frame(participant=rep(c(1:50), each = 400), trial=rep(c(1:200), each = 100), alternative = sample(1:4, 20000, T), attribute = sample(c("a","b","c","d"), 20000, T))

SSI = computeSSI(df, df, "participant", "trial", "alternative", "attribute", 4, 4, 4, 1000)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
