\name{RPB}
\alias{RPB}
\title{Robust Parametric Bootstrap Test based on MML estimators}
\description{
Computes the p-value of the robust parametric bootstrap \code{(RPB)} test for the equality of means of several long-tailed symmetric (LTS) distributions when the variances are unknown and arbitrary.}
\usage{
RPB(formula, data, alpha , verbose = TRUE, p_shape, repn)
}
\arguments{
  \item{formula}{a formula of the form left-hand-side\code{(lhs)} ~ right-hand-side\code{(rhs)}. \code{lhs} shows the observed values and \code{rhs} shows the group corresponding to the observed values.}
  \item{data}{data frame containing the variables in the formula.}
  \item{alpha}{the level of significance. Default is set to alpha = 0.05.}
  \item{verbose}{a logical for printing output to R console.}
  \item{p_shape}{shape parameter of the LTS distribution.}
  \item{repn}{replication number for performing the \code{RPB} test.}
}

\details{
\code{RPB} test based on modifed maximum likelihood (MML) estimators is proposed as a robust alternative to parametric bootstrap (PB) test proposed by Krishnamoorthy et. al (2007). See also Tiku (1967, 1968) for the details of MML estimators. The p-value for the \code{RPB} test is based on the replication number in the algorithm given by Guven et. al (2022).}

\value{
A list with class "\code{htest}" containing the following components:
  \item{p.value}{the p-value for the \code{RPB} test.}
  \item{alpha}{the level of significance.}
  \item{method}{a character string "Robust Parametric Bootstrap Test based on MML Estimators" indicating which test is used.}
  \item{data}{a data frame containing the variables.}
  \item{formula}{a formula of the form left-hand-side\code{(lhs)} ~ right-hand-side\code{(rhs)}. \code{lhs} shows the observed values and \code{rhs} shows the group corresponding to the observed values.}
}

\references{
G. Guven, S. Acitas, and B. Senoglu, B. RobustANOVA: An R Package for one-way ANOVA under heteroscedasticity and nonnormality. \emph{Under review}, 2022.

K. Krishnamoorthy, F. Lu, and T. Mathew. A parametric bootstrap approach for anova with unequal
variances: Fixed and random models. \emph{Computational Statistics & Data Analysis}, 51(12): 5731-5742,2007.

M. L. Tiku. Estimating the mean and standard deviation from a censored normal sample. \emph{Biometrika}, 54:155-165, 1967.

M. L. Tiku. Estimating the parameters of log-normal distribution from censored samples. \emph{Journal of the American Statistical Association}, 63(321): 134-140, 1968.

}

\author{Gamze Guven <gamzeguven@ogu.edu.tr>}


\examples{
library(RobustANOVA)

RPB(obs ~ methods, data = peak_discharge, alpha = 0.05, verbose = TRUE, p_shape=2.3, repn=5000)



}
\keyword{file}
