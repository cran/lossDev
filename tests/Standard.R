##set the working directory so the charts will be output into the correct directory

library(lossDev)
lossDevOptions(keepCodaOnDisk=FALSE)

## load the triangle, "data" loads the triangle as a data.frame so it must be coerced into a matrix
data(IncrementalGeneralLiablityTriangle)
IncrementalGeneralLiablityTriangle <- as.matrix(IncrementalGeneralLiablityTriangle)
IncrementalGeneralLiablityTriangle[1,1] <- NA


## load the stochastic inflation series, "data" loads the series as a data.frame so it must be coerced into a vector
data(CPI)
CPI <- as.matrix(CPI)[,1]
CPI.rate <- CPI[-1] / CPI[-length(CPI)] - 1

##restrict the cpi to only those years available when the triangle created
CPI.rate <- CPI.rate[as.integer(names(CPI.rate)) <= max(as.integer(dimnames(IncrementalGeneralLiablityTriangle)[[1]]))]


mi <- makeStandardAnnualInput(incremental.payments = IncrementalGeneralLiablityTriangle,
                              stoch.inflation.weight = 1,
                              non.stoch.inflation.weight = 0,
                              stoch.inflation.rate = CPI.rate,
                              exp.year.type = 'ay',
                              extra.dev.years=5,
                              use.ar1.in.calendar.year=FALSE,
                              use.ar1.in.exposure.growth=FALSE,
                              use.skew.t=TRUE)



mo <- runLossDevModel(mi,
                      burnIn=1000,
                      sampleSize=1000,
                      thin=1)


##check fit
predictedPayments(mo, plot=FALSE)
finalCumulativeDiff(mo, plot=FALSE)

triResi(mo,timeAxis='dy', plot=FALSE)
triResi(mo,timeAxis='cy', plot=FALSE)
triResi(mo,timeAxis='ey', plot=FALSE)

##development time
consumptionPath(mo, plot=FALSE)
consumptionPathTracePlot(mo, plot=FALSE)
numberOfKnots(mo, plot=FALSE)
rateOfDecay(mo, plot=FALSE)

standardDeviationVsDevelopmentTime(mo, plot=FALSE)


##exposure time
exposureGrowth(mo, plot=FALSE)
##exposureGrowthTracePlot(mo, plot=FALSE)
##meanExposureGrowth(mo, plot=FALSE)

##calendar time
calendarYearEffect(mo, restrictedSize=TRUE, plot=FALSE)
calendarYearEffectErrors(mo, plot=FALSE)
##calendarYearEffectErrorTracePlot(mo, plot=FALSE)

##inflation
##stochasticInflationRhoParameter(mo, plot=FALSE)
##stochasticInflationStationaryMean(mo, plot=FALSE)
##stochasticInflation(mo, plot=FALSE)


##extra
##skewnessParameter(mo, plot=FALSE)
##degreesOfFreedom(mo, plot=FALSE)
##scaleParameter(mo, plot=FALSE)



