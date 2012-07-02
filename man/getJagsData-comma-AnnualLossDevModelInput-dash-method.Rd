\name{getJagsData,AnnualLossDevModelInput-method}
\alias{getJagsData,AnnualLossDevModelInput-method}
\title{A method to create all the needed JAGS input common to both the standard model and the break model.}
\description{A method to create all the needed JAGS input common to both the standard model and the break model. Intended for internal use only.}
\details{There are currenly two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many data elements in common and the code to create these common elements is placed in this method.
The derived classes \code{StandardAnnualAggLossDevModelInput} and \code{BreakAnnualAggLossDevModelInput} should call this method via \code{NextMethod()} and then should append any needed data.

The returned value is a list containing the following elements with the following meanings:
\describe{
\item{\code{log.inc}}{A square matrix with only the upper right containing non-missing values of the log of the incremental payments.}
\item{\code{log.inc.index}}{A 2-column matrix.  Each row corresponds to a non-missing value in \code{log.inc}.  The first column gives the row of the non-missing value, the second the column.}
\item{\code{log.inc.index.length}}{The number of rows in \code{log.inc.index}.}
\item{\code{smooth.tau.h.2.log}}{A 2-valued vector giving the smoothing parameters for the changing (in development time) variance and skewness.}
\item{\code{h.same.as.last}}{A vector of 1's and 0's equal in length to the columns in \code{log.inc}.  A value of one means the variance and skewness in the corresponding column should be the same as in the previous column.}
\item{\code{K}}{Dimension of \code{log.inc}.}
\item{\code{H}}{Number of additional rows to estimate.}
\item{\code{L.vec}}{Vector equal in length to the number of total rows (observed plus forecast).  Each value is the column number of  the final non-zero incremental payment.}
\item{\code{delta.p}}{Vector equal in length to the number of total rows (observed plus forecast).  First column is (on the log scale) the mean of the final decay rate.  Second is standard deviation.}
\item{\code{v}}{Matrix with number of rows equal to total rows (observed plus forecast) and number of columns equal to total development years (observed plus forecast) minus the number of observed development years.}
\item{\code{ar1}}{Zero if calendar year effect does not have an ar1 component; 1 if it does.}
\item{\code{rho.prior}}{The parameters for the beta prior of the autoregressive parameter in calendar year effect.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{ar1.eta}}{Zero if exposure growth does not have an ar1 component; 1 if it does.}
\item{\code{rho.eta.prior}}{The parameters for the beta prior of the autoregressive parameter in exposure growth.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{allow.for.skew}}{Zero tells the model to assume zero skewness, one to use the skewed \eqn{t}.}
\item{\code{precision.for.skewness}}{The prior precision for the skewness parameter.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{df.for.skewness}}{The prior df for the skewness parameter.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{bounds.for.skewness}}{The the skewness parameter is bounted to prevent chains from getting stuck.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{precision.for.eta.mu}}{The prior precision for \code{eta.mu}.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{df.k}}{The parameter for the chi-sqr prior on the degrees of freedom.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{df.bounds}}{The lower and upper bounds on the degrees of freedom.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{a.ou.prior}}{The parameters for the beta prior of the autoregressive parameter for the stochastic inflation. Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{precision.for.b.ou}}{The prior precision for the intercept in the stochastic inflation ar1 process.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{sigma.eta.bounds}}{The lower and upper bounds on the standard deviation of the exposure year growth rate.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{sigma.kappa.bounds}}{The lower and upper bounds on the standard deviation of the calendar year growth rate. Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
\item{\code{stoch.log.inf}}{A vector giving the log inflation rate; missing future values are coded as NA.}
\item{\code{w.stoch.pct.inf.}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal to total development years. The portion of dollars inflating.}
\item{\code{non.stoch.log.inf}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal total development years. Non-stochastic inflation rate.  Cannot be \code{NA}.}
\item{\code{w.non.stoch.pct.inf}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal to total development years. The percent of dollars inflating.}
\item{\code{P}}{Single value.  The position in \code{stoch.log.inf.c} corresponding to the final observed diagonal in \code{log.inc}.}
\item{\code{stoch.log.inf.lower.bound}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal total development years. Floor for inflation rate.}
\item{\code{stoch.log.inf.upper.bound}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal total development years. Ceiling for inflation rate.}
\item{\code{estimate.a.ou}}{Single value (zero or one). Should the auto correlation coefficient for the Ornstein-Uhlenbeck process be estimated or should \code{fixed.a.ou} be used?}
\item{\code{fixed.a.ou}}{Single value. Possible non-stochastic  auto-correlation coefficient for the Ornstein-Uhlenbeck process.}
\item{\code{estimate.b.ou}}{Single value (zero or one). Should the stochastic rate of inflation have an estimated constant term or should \code{fixed.b.ou} be used?}
\item{\code{fixed.b.ou}}{Single value. Possible non-stochastic constant term for stochastic inflation rate.}
\item{\code{stoch.log.inf.known.mu}}{Single value.  Added to the log stochastic inflation rate after the ar1 estimation process.}
}}
\value{A named list of the specific model elements.  See Details for more information.}
\docType{methods}
\arguments{\item{object}{An object of type \code{AnnualAggLossDevModelInput} from which to collect the needed model input.}}
