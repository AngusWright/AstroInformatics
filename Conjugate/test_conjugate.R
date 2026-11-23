app_dir <- if (file.exists("Conjugate/app.R")) "Conjugate" else "."
source(file.path(app_dir, "conjugate.R"))

equal <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance)))
}

normal <- conjugate_update("normal", 4, 3, 1, 2, sigma = 2)
equal(normal$posterior$mean, 2.6)
equal(normal$posterior$sd^2, 0.8)
equal(normal$prior_weight, 0.2)
equal(1 / normal$posterior$sd^2, 1 / 4 + 4 / 4)

binomial <- conjugate_update("binomial", 10, 7, 2, 5)
equal(binomial$posterior$mean, 9 / 17)
equal(binomial$posterior$quantile(c(0.025, 0.975)), qbeta(c(0.025, 0.975), 9, 8))
poisson <- conjugate_update("poisson", 10, 40, 3, 2)
equal(poisson$posterior$mean, 43 / 12)
equal(poisson$posterior$sd^2, 43 / 12^2)
gamma <- conjugate_update("gamma", 5, 10, 3, 2, shape = 3)
equal(gamma$posterior$mean, 18 / 12)
equal(gamma$posterior$sd^2, 18 / 12^2)

# Independently check Bayes' product against actual observation likelihoods.
cases <- list(
  list(model = "normal", data = c(1, 2, 3, 6), statistic = 3,
       log_likelihood = function(theta) sum(dnorm(c(1, 2, 3, 6), theta, 2, log = TRUE))),
  list(model = "binomial", data = c(1, 0, 1, 1), statistic = 3,
       log_likelihood = function(theta) dbinom(3, 4, theta, log = TRUE)),
  list(model = "poisson", data = c(0, 2, 3, 1), statistic = 6,
       log_likelihood = function(theta) sum(dpois(c(0, 2, 3, 1), theta, log = TRUE))),
  list(model = "gamma", data = c(0.5, 1, 2, 2.5), statistic = 6,
       log_likelihood = function(theta) sum(dgamma(c(0.5, 1, 2, 2.5), 3, rate = theta, log = TRUE)))
)
for (case in cases) {
  result <- conjugate_update(case$model, 4, case$statistic, 2, 3, sigma = 2, shape = 3)
  x <- if (case$model == "binomial") c(0.2, 0.4, 0.7) else c(0.5, 1, 2)
  log_product <- log(result$prior$density(x)) + vapply(x, case$log_likelihood, numeric(1))
  log_posterior <- log(result$posterior$density(x))
  equal(diff(log_product), diff(log_posterior))
  equal(result$posterior$mean,
        result$prior_weight * result$prior$mean + (1 - result$prior_weight) * result$estimate)
  bounds <- if (case$model == "normal") c(-Inf, Inf) else if (case$model == "binomial") c(0, 1) else c(0, Inf)
  for (distribution in list(result$prior, result$likelihood, result$posterior)) {
    equal(integrate(distribution$density, bounds[1], bounds[2])$value, 1, tolerance = 1e-6)
  }
  empty <- conjugate_update(case$model, 0, 0, 2, 3, sigma = 2, shape = 3)
  equal(empty$posterior$density(x), empty$prior$density(x))
  stopifnot(is.null(empty$likelihood), is.na(empty$estimate), empty$prior_weight == 1)
}

for (successes in c(0, 10000)) {
  result <- conjugate_update("binomial", 10000, successes, 0.5, 0.5)
  stopifnot(all(is.finite(result$posterior$quantile(c(0.001, 0.999)))),
            result$posterior$mean > 0, result$posterior$mean < 1)
}
zero_counts <- conjugate_update("poisson", 10000, 0, 0.5, 0.1)
stopifnot(is.finite(zero_counts$posterior$mean), zero_counts$posterior$mean > 0)

large_n <- conjugate_update("normal", 1e7, 1, -1, 1, sigma = 2)
weak_prior <- conjugate_update("normal", 10, 1, -1, 1e6, sigma = 2)
equal(large_n$posterior$mean, 1, tolerance = 1e-6)
equal(large_n$posterior$sd^2, 4 / 1e7, tolerance = 1e-6)
equal(weak_prior$posterior$mean, 1)
equal(weak_prior$posterior$sd^2, 4 / 10)

bad_inputs <- list(
  list("normal", 10, 1, 0, 0),
  list("normal", 10, 1, 0, 1, sigma = 0),
  list("binomial", 10, 11, 1, 1),
  list("binomial", 10, 0, 0, 1),
  list("poisson", 10, 0.5, 1, 1),
  list("poisson", 0, 1, 1, 1),
  list("gamma", 10, 0, 1, 1),
  list("gamma", 10, 1, 1, 1, shape = 0)
)
for (args in bad_inputs) {
  stopifnot(inherits(try(do.call(conjugate_update, args), silent = TRUE), "try-error"))
}
cat("Analytic updates, Bayes products, normalisation, limits, and edge cases passed.\n")

app <- new.env()
source(file.path(app_dir, "app.R"), local = app, chdir = TRUE)
shiny::testServer(app$server, {
  session$setInputs(
    model = "normal", no_data = FALSE, log_n = 1, sample_mean = 1,
    sigma = 2, prior_mean = -1, prior_sd = 1, success_fraction = 0.7,
    beta_a = 2, beta_b = 5, count_mean = 4, gamma_mean = 2,
    shape = 3, gamma_a = 3, gamma_b = 2
  )
  for (model in c("normal", "binomial", "poisson", "gamma")) {
    session$setInputs(model = model, no_data = FALSE, log_n = 1)
    original <- analysis()
    session$setInputs(prior_mean = input$prior_mean + 0.5,
                      beta_a = input$beta_a + 1, gamma_a = input$gamma_a + 1)
    equal(analysis()$statistic, original$statistic)
    equal(analysis()$likelihood$density(c(0.25, 0.5, 0.75)),
          original$likelihood$density(c(0.25, 0.5, 0.75)))
    stopifnot(analysis()$posterior$mean != original$posterior$mean)
    stopifnot(
      nchar(output$summary) > 0, nchar(output$sufficiency) > 0,
      nchar(output$model_description$html) > 0,
      nchar(output$update_rule$html) > 0,
      nchar(output$information$html) > 0,
      nchar(output$role_explanation$html) > 0,
      nchar(output$density_plot$src) > 0,
      nchar(output$role_plot$src) > 0
    )
    session$setInputs(no_data = TRUE)
    equal(analysis()$posterior$mean, analysis()$prior$mean)
    stopifnot(nchar(output$density_plot$src) > 0,
              grepl("No observations", output$information$html))
  }
  for (fraction in c(0, 1)) {
    session$setInputs(model = "binomial", no_data = FALSE, log_n = 4,
                      success_fraction = fraction, beta_a = 0.5, beta_b = 0.5)
    stopifnot(nchar(output$density_plot$src) > 0, nchar(output$role_plot$src) > 0)
  }
  session$setInputs(model = "poisson", count_mean = 0, gamma_a = 0.5, gamma_b = 0.1)
  stopifnot(nchar(output$density_plot$src) > 0, nchar(output$role_plot$src) > 0)
  session$setInputs(model = "gamma", gamma_mean = 0.1, shape = 10, gamma_b = 20)
  stopifnot(nchar(output$density_plot$src) > 0, nchar(output$role_plot$src) > 0)
  session$setInputs(model = "normal", sigma = 0.1, prior_sd = 0.05, prior_mean = -5, sample_mean = 5)
  stopifnot(nchar(output$density_plot$src) > 0, nchar(output$role_plot$src) > 0)
})
cat("All four Shiny models, plots, tables, no-data states, and boundary controls passed.\n")
