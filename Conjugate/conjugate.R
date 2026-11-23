normal_distribution <- function(mean, sd) {
  list(
    label = sprintf("Normal(mean = %.4g, SD = %.4g)", mean, sd),
    mean = mean, sd = sd,
    density = function(x) dnorm(x, mean, sd),
    quantile = function(p) qnorm(p, mean, sd)
  )
}

beta_distribution <- function(a, b) {
  list(
    label = sprintf("Beta(%.4g, %.4g)", a, b),
    mean = a / (a + b),
    sd = sqrt(a * b / ((a + b)^2 * (a + b + 1))),
    density = function(x) dbeta(x, a, b),
    quantile = function(p) qbeta(p, a, b)
  )
}

gamma_distribution <- function(a, b) {
  list(
    label = sprintf("Gamma(shape = %.4g, rate = %.4g)", a, b),
    mean = a / b, sd = sqrt(a) / b,
    density = function(x) dgamma(x, shape = a, rate = b),
    quantile = function(p) qgamma(p, shape = a, rate = b)
  )
}

conjugate_update <- function(model, n, statistic, prior1, prior2,
                             sigma = 1, shape = 1) {
  stopifnot(
    model %in% c("normal", "binomial", "poisson", "gamma"),
    length(n) == 1, is.finite(n), n >= 0, n == round(n),
    length(statistic) == 1, is.finite(statistic),
    length(prior1) == 1, is.finite(prior1),
    length(prior2) == 1, is.finite(prior2), prior2 > 0
  )
  if (model != "normal") stopifnot(prior1 > 0, statistic >= 0)
  if (model != "normal" && n == 0) stopifnot(statistic == 0)
  if (model == "binomial") stopifnot(statistic <= n, statistic == round(statistic))
  if (model == "poisson") stopifnot(statistic == round(statistic))
  if (model == "gamma") {
    stopifnot(length(shape) == 1, is.finite(shape), shape > 0,
              if (n > 0) statistic > 0 else statistic == 0)
  }
  if (model == "normal") {
    stopifnot(length(sigma) == 1, is.finite(sigma), sigma > 0)
    prior <- normal_distribution(prior1, prior2)
    d <- 1 / prior2^2
    precision <- d + n / sigma^2
    weight <- d / precision
    posterior <- normal_distribution(
      weight * prior1 + (1 - weight) * statistic, sqrt(1 / precision)
    )
    likelihood <- if (n > 0) normal_distribution(statistic, sigma / sqrt(n)) else NULL
    estimate <- if (n > 0) statistic else NA_real_
    description <- sprintf("n = %d; observed sample mean = %.4g; known sigma = %.4g.", n, statistic, sigma)
    weight_description <- "Normal mean: prior weight = d/(d + n/sigma^2), where d = 1/s0^2."
  } else if (model == "binomial") {
    prior <- beta_distribution(prior1, prior2)
    posterior <- beta_distribution(prior1 + statistic, prior2 + n - statistic)
    likelihood <- if (n > 0) beta_distribution(statistic + 1, n - statistic + 1) else NULL
    weight <- (prior1 + prior2) / (prior1 + prior2 + n)
    estimate <- if (n > 0) statistic / n else NA_real_
    description <- sprintf("n = %d trials; s = %d successes; n - s = %d failures.", n, statistic, n - statistic)
    weight_description <- "Binomial probability: prior weight = (alpha + beta)/(alpha + beta + n). The prior concentration alpha + beta acts as an effective sample size for this mean update."
  } else if (model == "poisson") {
    prior <- gamma_distribution(prior1, prior2)
    posterior <- gamma_distribution(prior1 + statistic, prior2 + n)
    likelihood <- if (n > 0) gamma_distribution(statistic + 1, n) else NULL
    weight <- prior2 / (prior2 + n)
    estimate <- if (n > 0) statistic / n else NA_real_
    description <- sprintf("n = %d unit-exposure observations; total count = %d.", n, statistic)
    weight_description <- "Poisson rate: prior weight = beta/(beta + n). The prior rate beta acts as prior exposure in units of one observation."
  } else {
    prior <- gamma_distribution(prior1, prior2)
    posterior <- gamma_distribution(prior1 + n * shape, prior2 + statistic)
    likelihood <- if (n > 0) gamma_distribution(n * shape + 1, statistic) else NULL
    weight <- prior2 / (prior2 + statistic)
    estimate <- if (n > 0) n * shape / statistic else NA_real_
    description <- sprintf("n = %d; sum(x) = %.4g; known observation shape k = %.4g.", n, statistic, shape)
    weight_description <- "Gamma rate: prior weight = beta/(beta + sum(x)). The data estimate of rate is n*k/sum(x), not the sample mean."
  }
  if (n == 0) description <- "n = 0; no observed data. The posterior equals the prior."
  list(
    n = n, statistic = statistic, prior = prior, posterior = posterior,
    likelihood = likelihood, estimate = estimate, prior_weight = weight,
    data_description = description, weight_description = weight_description
  )
}
