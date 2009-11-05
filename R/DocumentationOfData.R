##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, National Council On Compensation Insurance Inc.,                        ##
##                                                                                              ##
##    This file is part of lossDev.                                                             ##
##                                                                                              ##
##    lossDev is free software: you can redistribute it and/or modify                           ##
##    it under the terms of the GNU General Public License as published by                      ##
##    the Free Software Foundation, either version 3 of the License, or                         ##
##    (at your option) any later version.                                                       ##
##                                                                                              ##
##    This program is distributed in the hope that it will be useful,                           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of                            ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             ##
##    GNU General Public License for more details.                                              ##
##                                                                                              ##
##    You should have received a copy of the GNU General Public License                         ##
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.                     ##
##                                                                                              ##
##################################################################################################

##' An and example of a cumulative loss triangle.
##'
##' Rows pertain to accident years; and columns pertain to development years.
##' For this triangle accident years range from 1974 to 1991 and development years range from 1 to 18.
##'
##' @references Hayne, Roger M. (2003), "Measurement of Reserve Variability," \emph{Casualty Actuarial Society Forum}, Fall 2003, pp. 141-172,
##' http://www.casact.org/pubs/forum/03fforum/03ff141.pdf.
##' @name CumulativeAutoBodilyInjuryTriangle
##' @usage CumulativeAutoBodilyInjuryTriangle
##' @docType data
NULL


##' Medical Care price index for the United States.
##'
##' @references \emph{bls.gov} \cr
##' \tabular{ll}{
##'   \bold{Not Seasonally Adjusted} \tab \cr
##'   \bold{Series Id:}              \tab CUUR0000SAM \cr
##'   \bold{Area:}                   \tab U.S. city average \cr
##'   \bold{Item:}                   \tab Medical care \cr
##'   \bold{Base Period:}            \tab 1982-84=100 \cr
##' }
##'
##' @name MCPI
##' @usage MCPI
##' @docType data
NULL

##' An and example of an incremental incurred loss triangle.
##'
##' Rows pertain to accident years; and columns pertain to development years.
##' For this triangle accident years range from 1981 to 1990 and development years range from 1 to 10.
##'
##' The triangle contains incurred incrementals of Automatic Facultative business in General Liability (excluding Asbestos & Environmental).
##'
##' @references
##' \emph{Historical Loss Development Study}, 1991 edition, published by the Reinsurance Association of America (RAA),
##' page 96. Quoted from, Thomas Mack (1995),
##' "Which Stochastic Model is Underlying the Chain Ladder Method," \emph{Casualty Actuarial Society Forum}, Fall 1995, pp. 229-240,
##' http://www.casact.org/pubs/forum/95fforum/95ff229.pdf.
##' @name IncrementalGeneralLiablityTriangle
##' @usage IncrementalGeneralLiablityTriangle
##' @docType data
NULL


##' Consumer price index (CPI-U) for the United States.
##'
##' @references \emph{bls.gov} \cr
##' \tabular{ll}{
##'   \bold{Not Seasonally Adjusted} \tab \cr
##'   \bold{Series Id:}              \tab CUUR0000SA0 \cr
##'   \bold{Area:}                   \tab U.S. city average \cr
##'   \bold{Item:}                   \tab All items \cr
##'   \bold{Base Period:}            \tab 1982-84=100 \cr
##' }
##'
##' @name CPI
##' @usage CPI
##' @docType data
NULL
