\name{RW}
\alias{RW}
\title{Robust Welch Test based on MML Estimators}
\description{
Computes the observed value of robust Welch \code{(RW)} test, degrees of freedoms (numerator and denominator) and the corresponding p-value for the equality of means of several long-tailed symmetric (LTS) distributions when the variances are unknown and arbitrary.}
\usage{
RW(formula, data, alpha=0.05, verbose = TRUE, p_shape)
}
\arguments{
  \item{formula}{a formula of the form left-hand-side\code{(lhs)} ~ right-hand-side\code{(rhs)}. \code{lhs} shows the observed values and \code{rhs} shows the group corresponding to the observed values.}
  \item{data}{data frame containing the variables in the formula.}
  \item{alpha}{the level of significance. Default is set to alpha = 0.05.}
  \item{verbose}{a logical for printing output to R console.}
  \item{p_shape}{shape parameter of the LTS distribution}
}

\details{
\code{RW} test based on modifed maximum likelihood (MML) estimators is proposed as a robust alternative to Welch's F test (Welch, 1951). The test statistic is formulated as follows

\deqn{RW= \frac{T(\hat{\mu}_1, \dots, \hat{\mu}_a;\hat{\sigma}_1^{2},\dots,\hat{\sigma}_a^{2})/(a-1)}{1+(2(a-2)/(3\nu_1))} }
where
\deqn{T(\hat{\mu}_1,\dots,\hat{\mu}_a; \hat{\sigma}_1^{2},\dots,\hat{\sigma}_a^{2})=\sum\limits_{i=1}^a \frac{M_i}{\hat{\sigma}_i^{2}} \hat{\mu}_i^{2}- \frac{(\sum\limits_{i=1}^a M_i\hat{\mu}_i/\hat{\sigma}_i^{2})^2}{\sum\limits_{i=1}^a M_i/\hat{\sigma}_i^{2}},}
\deqn{\nu_1= [\frac{3}{a^2-1} \sum\limits_{i = 1}^a \frac{1}{n_i-1}(1-( M_i/\hat{\sigma}_i^2)/( \sum\limits_{j= 1}^a  M_j/\hat{\sigma}_j^2))^{2}]^{-1},}
 \eqn{\hat{\mu}_{i}} and \eqn{\hat{\sigma}_{i}} (i=1,2,...,a) are the MML estimators of the location and scale parameters, respectively, see Tiku (1967, 1968) for the details of MML estimators.

The null hypothesis is rejected if the computed \code{RW} statistic is higher than the \eqn{(1-\alpha)}th quantile of the F distribution with a-1 and \eqn{\nu_{1}} degrees of freedom.

For further details, see Guven et al. (2022).

}

\value{
A list with class "\code{htest}" containing the following components:
  \item{statistic}{the observed value of the \code{RW} test statistic.}
  \item{dfs}{the numerator and the denominator degrees of freedom of the approximate F distribution.}
  \item{p.value}{the p-value for the \code{RW} test.}
  \item{alpha}{the level of significance.}
  \item{method}{a character string "Robust Welch Test based on MML Estimators" indicating which test is used.}
  \item{data}{a data frame containing the variables.}
  \item{formula}{a formula of the form left-hand-side\code{(lhs)} ~ right-hand-side\code{(rhs)}. \code{lhs} shows the observed values and \code{rhs} shows the group corresponding to the observed values.}
}

\references{
G. Guven, S. Acitas, and B. Senoglu, B. RobustANOVA: An R Package for one-way ANOVA under heteroscedasticity and nonnormality. \emph{Journal of Statistical Software}, 2022 (submitted).

M. L. Tiku. Estimating the mean and standard deviation from a censored normal sample. \emph{Biometrika}, 54:155-165, 1967.

M. L. Tiku. Estimating the parameters of log-normal distribution from censored samples. \emph{Journal of the American Statistical Association}, 63(321): 134-140, 1968.

B. L. Welch. On the comparison of several mean values: an alternative approach. \emph{Biometrika}, 38(3): 330-336, 1951.


}

\author{Gamze Guven <gamzeguven@ogu.edu.tr>}


\examples{
library(RobustANOVA)

RW(obs ~ methods, data = peak_discharge, alpha = 0.05, verbose = TRUE, p_shape=2.3)



}
\keyword{file}
