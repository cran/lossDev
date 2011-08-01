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
NULL




##' A class to hold \acronym{JAGS} output.  This class is only used internally.  No user-level function should return a class of this type.
##'
##' \code{NodeOutput} is a wrapper class for \code{mcarray}.
##' It is used to provide easy access to summary statistics.
##' Current slots are:
##' \describe{
##'   \item{\code{get.value.end}}{
##'     An environment containing a parameterless function called \code{get.value} which when called will return the mcarray for the node.
##'     It also contains \code{value.name} which is the name of the key (or file on the disk) if the value is stored on the disk.
##'   }
##'   \item{\code{mean}}{
##'     An array that is the marginalized mean of the value returned by calling \code{get.value}.
##'   }
##'   \item{\code{median}}{
##'     An array that is the marginalized median of the value returned by calling \code{get.value}.
##'   }
##'   \item{\code{sd}}{
##'     An array that is the marginalized standard deviation of the value returned by calling \code{get.value}.
##'   }
##' }
##' @name NodeOutput-class
##' @docType class
##' @seealso \code{\link{newNodeOutput}}.
setClass(
         'NodeOutput',
         representation(
                        get.value.env='environment',
                        mean='array',
                        median='array',
                        sd='array'
                        ))





##' A method to construct new object of type \code{NodeOutput}. Intended for internal use only.
##'
##' This method will return a valid NodeOutput object.
##' @param mcarray An S3 object of type \code{mcarray}.
##' @return An object of class \code{NodeOutput}.
##' @seealso \code{\linkS4class{NodeOutput}}
##  #import filehash only do this in zzz.R
newNodeOutput <- function(mcarray)
{
    ans <- new('NodeOutput')
    force(mcarray)

    ans@mean <- as.array(summary(mcarray, mean)[[1]])
    ans@median <- as.array(summary(mcarray, median)[[1]])
    my.sd <- function(x) sd(as.vector(x))[[1]]
    ans@sd <- as.array(summary(mcarray, my.sd)[[1]])



    ##create a new name for this object
    mutableState$CounterForCreatedCodas <- mutableState$CounterForCreatedCodas + 1
    value.name <- paste('object', mutableState$CounterForCreatedCodas, sep='')


    save.to.disk <- lossDevOptions()[['keepCodaOnDisk']]
    if(save.to.disk)
    {
        ans@get.value.env <-  new.env()
        ans@get.value.env$value.name <- value.name
        ans@get.value.env$get.value   <- function()
        {
            return(dbFetch(mutableState$fileHashDBForCodas, value.name))
        }

        reg.finalizer(ans@get.value.env,
                      function(env){

                          ##message('running reg.finalizer')
                          ##message(paste('value.name is', value.name))
                          ##message(paste('env$value.name is', env$value.name))
                          ##dbDelete(mutableState$fileHashDBForCodas, env$value.name)
                          dbDelete(mutableState$fileHashDBForCodas, value.name)
                      },
                      onexit=TRUE)

        dbInsert(mutableState$fileHashDBForCodas, value.name, mcarray)

        ##we have to remove this object from this closure, otherwise, since R must keep the closure around so it can still access "value.name"
        ##(or any other variable it doesn't know about for that matter)
        ##it will also keep "mcarray" around as well
        rm(mcarray)
    } else {

        ans@get.value.env <-  new.env()
        ans@get.value.env$get.value   <- function()
        {
            return(mcarray)
        }

    }



    #slot(ans, 'value') <- mcarray

    if(!validObject(ans))
        stop("could not create a valid")
    return(ans)
}


##' A method to override the behavoir of the function \code{slot}.
##'
##' In order to enhance the memory management, coda files are optionally stored on the harddrive in temporary files and loaded on an as needed basis.
##' By overriding this function, we are able to make this seamless.
##' Overriding the function \code{slot} is a slight abuse and, as such, may in the future be replaced by an accessor function.
##'
##' @name slot,NodeOutput,character-method
##' @param object The object of type \code{NodeOutput} with the slot to look up.
##' @param name A character value giving the name of the slot to look up.
##' @return Only if name is exactly \dQuote{value} will the method return the \code{mcarray} containing the coda.  Otherwise, it returns the result of \code{callNextMethod()}.
##' @docType methods
##' @seealso \code{\link{slot}}
setMethod('slot',
          signature(object='NodeOutput', name='character'),
          function(object, name)
      {
          if(length(name) != 1 || name != 'value')
              return(callNextMethod())

          return(object@get.value.env$get.value())


      })

