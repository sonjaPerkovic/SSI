\name{computeSSI}
\alias{computeSSI}
\title{A function that computes SSI
}
\description{
This function integrates all the necessary functions and computes the SSI.
It also creates a file with the SSI values for all the alternative- and
attribute-wise patterns found. The number of alternatives and attributes,
simulations as well as the minimum threshold for pattern length can be
adjusted by setting the function parameters.
}
\usage{
computeSSI(df, dfRan, participant, trial, alternative,
attribute, num_alt, num_att, threshold, iter)
}
\arguments{
  \item{df}{ object of class data frame
}
\item{dfRan}{ the same object which is used for creating random data
}
\item{participant}{ identifies each unique subject
}
\item{trial}{ identifies each unique trial
}
\item{alternative}{ column with eye fixations to different alternatives
}
\item{attribute}{ column with eye fixations to different attributes
}
\item{num_alt}{ number of alternatives in the experiment
}
\item{num_att}{ number of attributes in the experiment
}
\item{threshold}{ setting the threshold for pattern length to two or four
}
\item{iter}{ number of simulation iterations
}
}
\author{
Sonja Perkovic
}
\examples{
#IMPORTANT! Variables representing participants, trials, alternatives and attributes in your data should have the names\cr that match the ones in the example below!

dataSet <- data.frame(participant = rep(c(1:50), each = 400), trial = rep(c(1:200), each = 100), alternative = sample(1:4, 20000, TRUE), attribute = sample(c("a","b","c","d"), 20000, TRUE))

SSI <- computeSSI(dataSet, dataSet, "participant", "trial", "alternative", "attribute", 4, 4, 4, 10)
}
