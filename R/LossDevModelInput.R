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

##' @include zzz.R
NULL

##' The base input class for all models in \pkg{lossDev}.
##'
##' \code{LossDevModelInput} is the base input class for all model objects.
##' Derived classes must contain all needed data to construct input and initials for the \acronym{JAGS} model.
##' These are accessed via the S4 methods \dQuote{\code{\link{getJagsData}}} and \dQuote{\code{\link{getJagsInits}}.}
##' @name LossDevModelInput-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelOutput}}
setClass(
         'LossDevModelInput',
         representation(
                        #The model file to use.
                        modelFile='character',
                        #The name of the output class.
                        outputType='character',
                        'VIRTUAL'))


##' A generic function to run models in \pkg{lossDev}.
##'
##' Overriding methods must return a valid output object.
##' @name runLossDevModel
##' @param object The object containing the model to estimate.
##' @return object of class \code{LossDevModelOutput}.
##' @seealso \code{\link[=runLossDevModel,LossDevModelInput-method]{runLossDevModel("LossDevModelInput")}}
##' @exportMethod runLossDevModel
setGenericVerif('runLossDevModel',
                 function(object, ...)
                 standardGeneric('runLossDevModel'))




##' A method to run models in \pkg{lossDev}.
##'
##' This method returns a valid output object or flags an error.
##' This method is suitable for classes properly derived from class \code{LossDevModelInput} that properly override \dQuote{\code{\link{getJagsData}}} and \dQuote{\code{\link{getJagsInits}}}
##' and whose output type has a valid \code{getModelOutputNodes} method.
##'
##' \pkg{lossDev} sets the seed in each chain from a random number generated inside of \code{R}.
##' So to make a run reproducible, all one must do is set the seed (using \code{set.seed}) in \code{R} prior to the execution of this method.
##'
##' @name runLossDevModel,LossDevModelInput-method
##' @param object The object of type \code{LossDevModelInput} containing the model to estimate.
##' @param burnIn An integer to represent the number of initial \acronym{MCMC} iterations to be discarded. (The adaptive phase (\code{nAddapt}) is not considered part of \code{burnIn}.)
##' @param sampleSize An integer to represent the number of \acronym{MCMC} iterations to execute following the \code{burnIn}. (Actual number of iterations kept approximately \code{sampleSize / thin}.)
##' @param thin Keep every \code{thin}'th value of \code{sampleSize}.
##' @param nChains The number of \acronym{MCMC} chains to run.
##' @param nAddapt The length of the adaptive phase for the \acronym{MCMC} algorithm. (Default is \code{trunc(burnIn/4)+1}.)
##' @return An object of class \code{LossDevModelOutput}.
##' @docType methods
##' @seealso \code{\link{runLossDevModel}} \code{\link{set.seed}}
##  #import rjags only do this in zzz.R
setMethod(
          f='runLossDevModel',
          signature(object='LossDevModelInput'),
          definition=function(object, burnIn, sampleSize, thin=1, nChains=3, nAddapt=trunc(burnIn/4)+1)
      {
          time.begin<-Sys.time()   #call the time of your computer

          if(as.integer(burnIn) != burnIn || burnIn <= 0 || length(burnIn) != 1 )
              stop('Parameter "burnIn" must be a positive integer')
          if(as.integer(sampleSize) != sampleSize || sampleSize <= 0 || length(sampleSize) != 1)
              stop('Parameter "sampleSize" must be a positive integer')
          if(as.integer(thin) != thin || thin <= 0 || length(thin) != 1)
              stop('Parameter "thin" must be a positive integer')
          if(as.integer(nChains) != nChains || nChains <= 0 || length(nChains) != 1)
              stop('Parameter "nChains" must be a positive integer')
          if(as.integer(nAddapt) != nAddapt || nAddapt <= 0 || length(nAddapt) != 1)
              stop('Parameter "nAddapt" must be a positive integer')

          ##ans <- new(getOutputType(object))
          ans <- new(object@outputType)

          parameters.to.save. <- getModelOutputNodes(ans)
          ##print(parameters.to.save.)


          ##We will NOT count the addaptive phase as part of the burnin.
          ##ans@burnIn <- as.integer(burnIn)
          ##ans@sampleSize <- as.integer(sampleSize)
          ##ans@thin <- as.integer(thin)
          ##ans@nChains <- as.integer(nChains)

          message(paste('Preparing Jags Model\nadapting for', nAddapt, 'iterations\n\n'))
          gc()

          rngs <- rep(paste('base', c('Wichmann-Hill', 'Marsaglia-Multicarry', 'Super-Duper', 'Mersenne-Twister'), sep='::'), length.out=nChains)

          if(nChains > 100) stop(paste('Why are you running so many chains?', nChains, 'is too many.  Let\'s keep it below 100, OK?'))
          gen.seed <- function() ceiling(runif(1, 1, 10000))

          rng.seeds <- gen.seed()
          for(i in seq(1, nChains)[-1])
          {
              prop.seed <- gen.seed()
              while(prop.seed %in%  rng.seeds)
                  prop.seed <- gen.seed()

              rng.seeds[i] <- prop.seed
          }

          master.inits.f <- getJagsInits(object)

          inits.f <- function(chain)
          {
              ans <- master.inits.f()
              ans[['.RNG.name']] <- rngs[chain]
              ans[['.RNG.seed']] <- rng.seeds[chain]

              return(ans)

          }

          jm <- jags.model(file=file.path(myLibPath(), myPkgName(), 'models', object@modelFile),
                           data=getJagsData(object),
                           inits=inits.f,
                           n.chains=nChains,
                           n.adapt=nAddapt)


          message(paste('Burning-In Jags Model for', burnIn, 'iterations\n', 'Total Burn-In = ', burnIn))
          update(jm, burnIn)

          message(paste('Sampling Jags Model for', sampleSize, 'iterations Thin =', thin,'\n', 'This will result in ~', sampleSize / thin, 'Samples'))
          output <- jags.samples(jm, parameters.to.save., sampleSize, thin)
          rm(jm)
          gc()
          gc()


          ans@input <- object
          for(i in parameters.to.save.)
              slot(ans,i) <- newNodeOutput(output[[i]])
          ##slot(ans,i) <- new('NodeOutput', value=new('safe.mcmc', value=output[[i]]))
          rm(output)

          if(!validObject(ans))
              stop('A valid output could not be created')

          print(paste('Update took', format(Sys.time() - time.begin)))

          ##attributes(ans)[['jags.model']] <- jm
          return(invisible(ans))
      })


##' A method to retrieve initial values for a \acronym{JAGS} model. Intended for internal use only.
##'
##' All classes derived from class \code{LossDevModelInput} should override this method.
##' This the overriding method should return a parameterless function that, when evaluated, returns a named list of initial values.
##' @name getJagsInits
##' @param object The object from which to retrieve initial values.
##' @return A function. Overriding methods should return a parameterless function that, when evaluated, returns a named list of initial values.
setGenericVerif('getJagsInits',
                 function(object)
                 standardGeneric('getJagsInits'))

##' A method to retrieve data for a \acronym{JAGS} model. Intended for internal use only.
##'
##' Any classes derived from class \code{LossDevModelInput} must provide a method to override this generic function.
##' This overriding method should return a named list of data values.
##' @name getJagsData
##' @param object The object from which to retrieve initial values.
##' @return This method returns a named list of data values.
setGenericVerif('getJagsData',
                 function(object)
                 standardGeneric('getJagsData'))
