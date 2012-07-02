library(rjags)

simple.model <-
  "model
{
     y ~ dnorm(mu, tau)
     mu ~ dnorm(0, 1)
     tau ~ dgamma(0.1, 0.1)
}
"

simple.model.file <- tempfile()
cat(simple.model, file=simple.model.file)

jags.data <- list(y=0)
jm <- jags.model(file=simple.model.file,
                 data=jags.data)

update(jm, 100)
