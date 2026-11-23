shared <- new.env()
sys.source("TeachingTools/teaching.R", envir = shared)
load_model <- function(directory) {
  model <- new.env(parent = shared)
  sys.source(file.path(directory, "model.R"), envir = model)
  model
}
equal <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(unname(actual), unname(expected),
                            tolerance = tolerance, check.attributes = FALSE)))
}
seeded <- shared$lesson_seeded

selection <- load_model("SelectionEffects")
equal(selection$selection_observed(c(0, 1), c(100, 100), 16), c(FALSE, TRUE))
catalogue <- seeded(42, selection$calculate(list(limit = 19), 1))
stopifnot(abs(cor(catalogue$distance, catalogue$luminosity)) < 0.06,
          mean(catalogue$luminosity[catalogue$detected]) > 0.2)
deeper <- selection$selection_observed(catalogue$luminosity, catalogue$distance, 22)
stopifnot(all(deeper[catalogue$detected]))

eddington <- load_model("EddingtonBias")
for (noise in c(0, 0.1, 0.6)) {
  equal(integrate(function(x) eddington$mass_density(x, 1.8, noise), 4, 17,
                  subdivisions = 500)$value, 1, 1e-6)
}
mass <- seeded(42, eddington$calculate(list(slope = 1.8, noise = 0.4, threshold = 10.5), 1))
stopifnot(abs(mean(mass$observed - mass$true)) < 0.02, mass$up > mass$down)
equal(sum(mass$observed >= 10.5) - sum(mass$true >= 10.5), mass$up - mass$down)

hierarchy <- load_model("HierarchicalModels")
x <- c(-1, 1, 2)
error <- c(0.3, 0.8, 2)
prior_covariance <- diag(0.7^2, 3) + matrix(25, 3, 3)
posterior_covariance <- solve(solve(prior_covariance) + diag(1 / error^2))
posterior <- hierarchy$normal_hierarchy(x, error, 0.7)
equal(posterior$mean, drop(posterior_covariance %*% (x / error^2)))
equal(posterior$variance, diag(posterior_covariance))
pooled <- hierarchy$normal_hierarchy(x, error, 0)
equal(pooled$mean, rep(pooled$population_mean, 3))
equal(pooled$variance, rep(pooled$population_sd^2, 3))
independent <- hierarchy$normal_hierarchy(x, error, 1e6)
equal(independent$mean, x)
equal(independent$variance, error^2)

evidence <- load_model("ModelEvidence")
alternative <- integrate(function(a) dnorm(1.1, a, 0.5) * dnorm(a, 0, 2), -Inf, Inf)$value
score <- evidence$line_evidence(1.1, 0.5, 2)
equal(exp(score$log_bf), alternative / dnorm(1.1, 0, 0.5))
equal(score$probability / (1 - score$probability), exp(score$log_bf))
parameters <- list(amplitude = 0.5, noise = 0.3, prior = 1)
narrow <- seeded(42, evidence$calculate(parameters, 1))
parameters$prior <- 5
broad <- seeded(42, evidence$calculate(parameters, 1))
equal(narrow$flux, broad$flux)
equal(narrow$evidence$log_maximum_ratio, broad$evidence$log_maximum_ratio)
stopifnot(evidence$line_evidence(0, 0.5, 5)$log_bf < evidence$line_evidence(0, 0.5, 1)$log_bf)

testing <- load_model("MultipleTesting")
equal(testing$search_false_probability(3, 0)$global, 0)
equal(testing$search_false_probability(3, 1)$global, testing$search_false_probability(3, 1)$single)
equal(testing$search_false_probability(3, 10, TRUE)$single, 2 * pnorm(-3))
false_alarm <- seeded(42, mean(colSums(matrix(rnorm(400000), 20) > 2.5) > 0))
stopifnot(abs(false_alarm - testing$search_false_probability(2.5, 20)$global) < 0.01)

marginal <- load_model("Marginalisation")
for (location in c(-1, 0, 1.3, 3)) {
  equal(integrate(function(y) marginal$mixture_joint(location, y, 0.7, 1.5), -Inf, Inf)$value,
        marginal$mixture_marginal(location, 0.7, 1.5), 1e-6)
}
equal(integrate(function(x) marginal$mixture_marginal(x, 0.7, 1.5), -Inf, Inf)$value, 1)
equal(integrate(function(x) marginal$mixture_marginal(x, 0.7, 1.5), -0.5, 0.8)$value,
      marginal$mixture_interval(-0.5, 0.8, 0.7, 1.5))

base_rates <- load_model("BaseRates")
ppv <- base_rates$positive_predictive_value(0.001, 0.95, 0.01)
equal(ppv, 950 / (950 + 9990))
stopifnot(ppv < 0.1,
          base_rates$positive_predictive_value(0.1, 0.95, 0.01) > ppv,
          base_rates$positive_predictive_value(0.001, 0.95, 0.001) > ppv)

prediction <- load_model("PosteriorPredictive")
small <- prediction$predictive_normal(c(1, 2, 3), 1, 2)
equal(small$mean, 1.875)
equal(small$variance, c(1.5625, 2.5625, 6.5625))
large <- prediction$predictive_normal(rep(2, 100000), 1, 2)
stopifnot(large$variance[1] < 0.0001)
equal(diff(large$variance), c(1, 4))

forward <- load_model("ForwardModelling")
encounters <- seeded(42, forward$prepare(list(impact = 2, speed = 1, uncertainty = 0)))
stopifnot(all(is.finite(encounters$nominal)))
for (path in encounters$ensemble) equal(path, encounters$nominal)
coarse <- forward$restricted_encounter(c(0.8, 1.2), c(0, 1), 1000, 1, steps = 80, dt = 0.05)
fine <- forward$restricted_encounter(c(0.8, 1.2), c(0, 1), 1000, 1, steps = 160, dt = 0.025)
stopifnot(max(abs(coarse[81, , ] - fine[161, , ])) < 0.01,
          max(abs(sqrt(rowSums(fine[161, , ]^2)) - c(0.8, 1.2))) < 0.002)

correlation <- load_model("Correlation")
galaxies <- seeded(42, correlation$calculate(list(scatter = 0.4, centre = 0.5, width = 0.1), 1))
stopifnot(abs(galaxies$full_correlation - 1 / (1 + 3 * 0.4^2)) < 0.05,
          abs(cor(galaxies$x + 2 * galaxies$distance, galaxies$y + 2 * galaxies$distance)) < 0.08)

survey <- load_model("SurveyShift")
original <- seeded(42, survey$calculate(list(offset = 0, noise = 0), 1))
shifted <- seeded(42, survey$calculate(list(offset = 1.5, noise = 0), 1))
equal(original$training, shifted$training)
equal(original$original_prediction, shifted$original_prediction)
equal(shifted$target - original$target, matrix(c(1.5, 0.6), 400, 2, byrow = TRUE))
stopifnot(mean(original$target_prediction == original$label) > 0.95,
          mean(shifted$target_prediction == shifted$label) < 0.8)
relabelled <- survey$fit_colour_classifier(original$training, 1 - original$label)
equal(survey$classify_colours(relabelled, original$original), 1 - original$original_prediction)

photons <- load_model("PhotonCounting")
zero <- seeded(42, photons$calculate(list(rate = 2.5, exposure = 4), 0))
short <- seeded(42, photons$calculate(list(rate = 2.5, exposure = 4), 0.5))
long <- seeded(42, photons$calculate(list(rate = 2.5, exposure = 4), 1))
stopifnot(!length(zero$arrivals), all(zero$counts == 0),
          all(diff(long$arrivals) > 0), all(long$arrivals <= long$exposure),
          all(long$counts >= short$counts),
          abs(mean(long$counts) - 10) < 0.25, abs(var(long$counts) - 10) < 0.7)
equal(short$arrivals, long$arrivals[long$arrivals <= short$exposure])

stacking <- load_model("CentralLimit")
for (population in c("exponential", "mixture", "cauchy")) {
  parameters <- list(population = population, repetitions = 5000, n = 100)
  experiments <- seeded(42, stacking$prepare(parameters))
  averages <- stacking$calculate(parameters, 1, experiments)
  if (population == "cauchy") {
    stopifnot(!averages$finite, is.na(averages$mean_sd))
    equal(quantile(averages$means, c(0.25, 0.5, 0.75), names = FALSE), c(-1, 0, 1), 0.15)
  } else {
    stopifnot(abs(mean(averages$statistic)) < 0.05,
              abs(var(averages$statistic) - 1) < 0.06)
    equal(averages$mean_sd, if (population == "mixture") sqrt(2.1875) / 10 else 0.1)
  }
}
window <- shared$lesson_histogram(c(-1, 0, 0.5, 1, 2, 3), c(0, 1, 2))
equal(window$density, c(1 / 3, 1 / 3))
equal(window$excluded, 2)
equal(shared$lesson_histogram(numeric(0), c(0, 1, 2))$density, c(0, 0))

inverse <- load_model("InverseSampling")
for (family in c("power", "exponential", "lognormal")) {
  for (slope in c(1, 2.35)) {
    distribution <- inverse$positive_distribution(family, slope)
    probability <- c(0.00001, 0.05, 0.5, 0.9, 0.99999)
    equal(distribution$cdf(distribution$quantile(probability)), probability)
    bounds <- if (family == "power") c(0.1, 20) else c(0, Inf)
    equal(integrate(distribution$density, bounds[1], bounds[2])$value, 1, 1e-6)
    if (family == "power") {
      equal(integrate(function(z) log(10) * 10^z * distribution$density(10^z),
                      log10(0.1), log10(20))$value, 1)
    }
  }
}

robust <- load_model("RobustStatistics")
moderate <- robust$robust_summaries(c(-2, -1, 0, 1, 100))
extreme <- robust$robust_summaries(c(-2, -1, 0, 1, 1000))
equal(moderate[c("median", "nmad")], extreme[c("median", "nmad")])
stopifnot(extreme["sd"] > 9 * moderate["sd"])
equal(extreme["mean"] - moderate["mean"], 180)

density_model <- load_model("DensityEstimation")
parameters <- list(n = 1000, origin = 0, width = 0.25, bandwidth = 0.15)
sharp <- seeded(42, density_model$calculate(parameters, 1))
parameters$origin <- 0.7
parameters$width <- 0.8
parameters$bandwidth <- 1
smooth <- seeded(42, density_model$calculate(parameters, 1))
equal(sharp$values, smooth$values)
equal(sum(smooth$histogram$density * diff(smooth$histogram$breaks)), 1)
equal(integrate(density_model$stellar_colour_density, -Inf, Inf)$value, 1)
stopifnot(length(sharp$peaks) >= 2, length(smooth$peaks) == 1)

cat("All 16 lessons: analytic probabilities, density normalisation, posterior covariance,\n",
    "simulation limits, sampling distributions and data-invariance checks passed.\n", sep = "")
