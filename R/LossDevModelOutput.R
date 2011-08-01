##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, 2010, 2011 National Council On Compensation Insurance Inc.,             ##
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

##' @include zzz.R
##' @include LossDevModelInput.R
NULL

##' The base output class for all models in \pkg{lossDev}.
##'
##' \code{LossDevModelOutput} is the base class for all output model objects.
##' Derived classes should contain all needed output from a \acronym{JAGS} run and the input object associated with that run in the slot \dQuote{input.}
##' Currenly, only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes must be the exact name of some settable node in the model.
##' This is because \code{getModelOutputNodes} currently looks at the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
##' Any class adding extra non-model node slots must also override the method \code{getModelOutputNodes}.
##' @name LossDevModelOutput-class
##' @docType class
##' @seealso \code{\linkS4class{AnnualAggLossDevModelInput}}
setClass(
         'LossDevModelOutput',
         representation(
                        input='LossDevModelInput',
                        'VIRTUAL'))

##' A method to determine which \acronym{JAGS} nodes the output object is expecting.  Intended for internal use only.
##'
##' This method is used to query the nodes a particular \code{LossDevModelOutput} is expecting.
##' @name getModelOutputNodes
##' @param object The object from which to get the expected nodes.
##' @return A character vector with each element corresponding to the name of an expected node.
##' @seealso \code{\link[=getModelOutputNodes,LossDevModelOutput-method]{getModelOutputNodes("LossDevModelOutput")}}
setGenericVerif('getModelOutputNodes',
                function(object)
                standardGeneric('getModelOutputNodes'))

##' A method to determine which \acronym{JAGS} nodes the output object is expecting.  Intended for internal use only.
##'
##' This method works by examining all the slot names in \code{object} and returning all of them except for a slot named \dQuote{input.}
##' @name getModelOutputNodes,LossDevModelOutput-method
##' @param object The object of type \code{LossDevModelOutput} from which to get the expected nodes.
##' @return A character vector with each element corresponding to the name of an expected node.
##' @docType methods
##' @seealso \code{\link{getModelOutputNodes}}
setMethod('getModelOutputNodes',
          signature(object='LossDevModelOutput'),
          function(object)
      {
          ans <- slotNames(class(object)[1])
          return(ans[-which(ans=='input')])
      })
