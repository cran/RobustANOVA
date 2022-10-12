\name{peak_discharge}
\alias{peak_discharge}
\title{Peak Disharge Data}
\usage{
peak_discharge
}
\description{
The "peak discharge data" first given by Montgomery (2005) consists of four different methods of estimating flood flow frequency.
}
\examples{
library(RobustANOVA)
peak_discharge$obs;
peak_discharge$methods;
}

\value{
  \item{obs}{Flood flow frequency (in cubic feet per second)}
  \item{methods}{Methods of estimating flood flow frequency.}
}
\author{Gamze Guven}
\references{
D. C. Montgomery. Design and analysis of experiments. John wiley & sons, 2005.
}

