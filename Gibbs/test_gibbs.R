app_dir <- if (file.exists("Gibbs/app.R")) "Gibbs" else "."
source(file.path(app_dir, "gibbs.R"))

equal <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance)))
}

gaussian <- gibbs_target("gaussian", location = c(1, -2), scale = c(2, 3), rho = 0.6)
conditional1 <- gaussian$conditional(c(2, -1), 1)
conditional2 <- gaussian$conditional(c(2, -1), 2)
equal(conditional1$quantile(0.5), 1.4)
equal(conditional1$density(1.4), dnorm(1.4, 1.4, 1.6))
equal(conditional2$quantile(0.5), -1.1)
equal(conditional2$density(-1.1), dnorm(-1.1, -1.1, 2.4))

student <- gibbs_target("student", location = c(1, -2), scale = c(2, 3), rho = 0.6, df = 5)
width <- 2 * sqrt((1 - 0.6^2) * (5 + 1 / 9) / 6)
equal(student$conditional(c(2, -1), 1)$density(1.4), dt(0, 6) / width)
poisson <- gibbs_target("poisson", n = 10, total = 30, a = 2, c = 2, d = 1)
equal(poisson$conditional(c(3, 2), 1)$quantile(0.5), qgamma(0.5, 32, rate = 12))
equal(poisson$conditional(c(3, 2), 2)$quantile(0.5), qgamma(0.5, 4, rate = 4))
mixture <- gibbs_target("mixture", location = c(1, -2), scale = c(2, 3),
                         rho = 0.6, separation = 3, weight = 0.3)
collapsed <- gibbs_target("mixture", location = c(1, -2), scale = c(2, 3),
                           rho = 0.6, separation = 0, weight = 0.2)
equal(collapsed$log_density(c(-1, 0, 2), c(2, 0, -1)),
      gaussian$log_density(c(-1, 0, 2), c(2, 0, -1)))
for (j in 1:2) {
  equal(collapsed$conditional(c(2, 1), j)$quantile(c(0.001, 0.5, 0.999)),
        gaussian$conditional(c(2, 1), j)$quantile(c(0.001, 0.5, 0.999)))
  density <- function(x) mixture$marginal(x, j)
  equal(integrate(density, -Inf, Inf)$value, 1, 1e-6)
  mean <- integrate(function(x) x * density(x), -Inf, Inf)$value
  variance <- integrate(function(x) (x - mean)^2 * density(x), -Inf, Inf)$value
  equal(mean, mixture$mean[j], 1e-6)
  equal(sqrt(variance), mixture$sd[j], 1e-6)
  conditional <- mixture$conditional(c(2, 1), j)
  for (p in c(0.001, 0.1, 0.5, 0.9, 0.999)) {
    equal(integrate(conditional$density, -Inf, conditional$quantile(p),
                    rel.tol = 1e-10, abs.tol = 1e-12)$value, p, 1e-6)
  }
  equal(conditional$quantile(c(0, 1)), c(-Inf, Inf))
}
separated <- gibbs_target("mixture", rho = 0, separation = 4)
thin_conditional <- gibbs_target("mixture", rho = -0.98, separation = 8)$conditional(c(0, 0), 1)
equal(thin_conditional$quantile(0.5), 0)
stopifnot(integrate(separated$conditional(c(0, -4), 1)$density, -Inf, 0)$value > 0.97,
          integrate(separated$conditional(c(0, 4), 1)$density, -Inf, 0)$value < 0.03,
          is.finite(mixture$log_density(1e6, -1e6)),
          all(is.finite(mixture$conditional(c(1e6, -1e6), 1)$quantile(c(0.001, 0.999)))))

# Every full conditional must agree with the corresponding slice of the joint.
for (model in list(gaussian, student, poisson, mixture)) {
  state <- c(2, 1)
  bounds <- if (model$positive) c(0, Inf) else c(-Inf, Inf)
  for (j in 1:2) {
    conditional <- model$conditional(state, j)
    equal(integrate(conditional$density, bounds[1], bounds[2])$value, 1, 1e-6)
    x <- c(0.5, 1, 2)
    joint <- if (j == 1) model$log_density(x, state[2]) else model$log_density(state[1], x)
    equal(diff(log(conditional$density(x))), diff(joint))
  }
}

# Verify the second draw uses the newly sampled first coordinate, not the old one.
set.seed(99)
expected <- matrix(NA_real_, 7, 2, dimnames = list(NULL, gaussian$axes))
expected[1, ] <- c(-2, 3)
for (sweep in 1:3) {
  before <- expected[2 * sweep - 1, ]
  x <- rnorm(1, 1 + 0.6 * 2 / 3 * (before[2] + 2), 1.6)
  y <- rnorm(1, -2 + 0.6 * 3 / 2 * (x - 1), 2.4)
  expected[2 * sweep, ] <- c(x, before[2])
  expected[2 * sweep + 1, ] <- c(x, y)
}
seed_before <- .Random.seed
chain <- gibbs_chain(gaussian, c(-2, 3), 3, 99)
equal(chain$states, expected)
stopifnot(identical(.Random.seed, seed_before))
equal(gibbs_chain(gaussian, c(-2, 3), 10, 99)$states[1:7, ], expected)
for (model in list(gaussian, student, poisson, mixture)) {
  generated <- gibbs_chain(model, c(2, 1), 20, 666)
  differences <- generated$states[-1, ] - generated$states[-nrow(generated$states), ]
  equal(differences[seq(1, 40, by = 2), 2], rep(0, 20))
  equal(differences[seq(2, 40, by = 2), 1], rep(0, 20))
  if (model$positive) stopifnot(all(generated$states > 0))
  first <- gibbs_frame(generated, 1)
  second <- gibbs_frame(generated, 2)
  equal(first$after, second$before)
  equal(nrow(gibbs_samples(generated, 1)), 0L)
  equal(nrow(gibbs_samples(generated, 2)), 1L)
  equal(gibbs_samples(generated, 5), generated$states[c(3, 5), ])
  equal(gibbs_samples(generated, 5, burn = 1), generated$states[5, , drop = FALSE])
  equal(nrow(gibbs_samples(generated, 5, burn = 3)), 0L)
}
cat("Exact conditionals, joint-density slices, sequential updates, RNG isolation, and sweep accounting passed.\n")

for (rho in c(-0.98, 0, 0.98)) {
  samples <- gibbs_chain(gibbs_target(rho = rho), c(-2.5, 2.5), 20, 1)$states
  stopifnot(all(is.finite(samples)))
  samples <- gibbs_chain(gibbs_target("mixture", rho = rho, separation = 8, weight = 0.05),
                         c(-2.5, 2.5), 20, 1)$states
  stopifnot(all(is.finite(samples)))
}
normal_samples <- gibbs_samples(
  gibbs_chain(gibbs_target(rho = 0.8), c(-2.5, 2.5), 30000, 666), 60000, 1000
)
stopifnot(max(abs(colMeans(normal_samples))) < 0.06,
          max(abs(apply(normal_samples, 2, sd) - 1)) < 0.05,
          abs(cor(normal_samples)[1, 2] - 0.8) < 0.025)
student_samples <- gibbs_samples(
  gibbs_chain(gibbs_target("student", rho = -0.5, df = 6), c(-2.5, 2.5), 30000, 123), 60000, 1000
)
stopifnot(max(abs(colMeans(student_samples))) < 0.08,
          max(abs(apply(student_samples, 2, sd) - sqrt(1.5))) < 0.1,
          abs(cor(student_samples)[1, 2] + 0.5) < 0.04)
poisson_samples <- gibbs_samples(gibbs_chain(poisson, c(1, 1), 20000, 42), 40000, 1000)
kernel <- function(x) exp(31 * log(x) - 10 * x - 4 * log1p(x))
normalisation <- integrate(kernel, 0, Inf)$value
expected_lambda <- integrate(function(x) x * kernel(x), 0, Inf)$value / normalisation
expected_beta <- integrate(function(x) 4 / (1 + x) * kernel(x), 0, Inf)$value / normalisation
stopifnot(abs(mean(poisson_samples[, 1]) - expected_lambda) < 0.04,
          abs(mean(poisson_samples[, 2]) - expected_beta) < 0.04)
mixing_target <- gibbs_target("mixture", location = c(0.4, -0.7), scale = c(0.8, 1.3),
                              rho = 0.2, separation = 2, weight = 0.35)
mixture_samples <- gibbs_samples(
  gibbs_chain(mixing_target, c(-2.5, 2.5), 20000, 42), 40000, 1000
)
stopifnot(max(abs(colMeans(mixture_samples) - mixing_target$mean)) < 0.07,
          max(abs(apply(mixture_samples, 2, sd) - mixing_target$sd)) < 0.08,
          abs(cov(mixture_samples)[1, 2] - 0.8 * 1.3 * (0.2 + 0.35 * 0.65 * 4)) < 0.1,
          abs(mean(mixture_samples[, 1] < 0.4) - (0.35 * pnorm(1) + 0.65 * pnorm(-1))) < 0.025)
zero_counts <- gibbs_chain(gibbs_target("poisson", total = 0, a = 1), c(1, 1), 30, 2)
stopifnot(all(is.finite(zero_counts$states)), all(zero_counts$states > 0))
stopifnot(
  inherits(try(gibbs_target(rho = 1), silent = TRUE), "try-error"),
  inherits(try(gibbs_target("student", df = 2), silent = TRUE), "try-error"),
  inherits(try(gibbs_target("mixture", separation = -1), silent = TRUE), "try-error"),
  inherits(try(gibbs_target("mixture", weight = 0), silent = TRUE), "try-error"),
  inherits(try(gibbs_target("mixture", weight = 1), silent = TRUE), "try-error"),
  inherits(try(gibbs_chain(poisson, c(0, 1)), silent = TRUE), "try-error")
)
cat("Stationary Gaussian/Student-t/mixture moments, quadrature-checked moments, and support boundaries passed.\n")

app <- new.env()
source(file.path(app_dir, "app.R"), local = app, chdir = TRUE)
shiny::testServer(app$server, {
  session$setInputs(
    family = "gaussian", location1 = 0, location2 = 0, scale1 = 1, scale2 = 1,
    rho = 0.8, df = 5, n = 10, total = 30, a = 2, hyper_shape = 2, hyper_rate = 1,
    separation = 4, weight = 0.5,
    sweeps = 20, seed = 666, start1 = -2.5, start2 = 2.5,
    start_lambda = 1, start_beta = 1, delay = "250", trail = 20, burn = 2,
    show_samples = TRUE, jump_to = 0
  )
  original <- chain()$states
  stopifnot(frame() == 0, !playing(), nrow(completed()) == 0,
            nchar(output$joint_plot$src) > 0, nchar(output$conditional_plot$src) > 0)
  session$setInputs(step = 1)
  stopifnot(frame() == 1, nrow(completed()) == 0, current()$coordinate == 1)
  session$setInputs(step = 2)
  stopifnot(frame() == 2, nrow(completed()) == 1, current()$coordinate == 2)
  session$setInputs(sweep = 1)
  stopifnot(frame() == 4, nrow(completed()) == 2)
  session$setInputs(rewind = 1)
  stopifnot(frame() == 0, !playing())

  session$setInputs(play = 1)
  stopifnot(playing(), frame() == 1)
  session$elapse(251)
  stopifnot(frame() == 2)
  session$setInputs(play = 2)
  paused <- frame()
  session$elapse(1000)
  stopifnot(!playing(), frame() == paused)
  session$setInputs(jump_to = 39, jump = 1)
  stopifnot(frame() == 39, !playing())
  session$setInputs(play = 3)
  stopifnot(frame() == 40, !playing())
  session$elapse(1000)
  stopifnot(frame() == 40)
  session$setInputs(play = 4)
  stopifnot(frame() == 1, playing())
  session$setInputs(step = 3)
  stopifnot(frame() == 2, !playing())
  session$setInputs(jump_to = 41, jump = 2)
  stopifnot(frame() == 2)

  session$setInputs(burn = 3, trail = 10, delay = "1000", show_samples = FALSE)
  equal(chain()$states, original)
  stopifnot(frame() == 2)
  session$setInputs(finish = 1)
  stopifnot(frame() == 40, nrow(retained()) == 17)
  for (family in c("gaussian", "student", "mixture", "poisson")) {
    session$setInputs(family = family)
    session$setInputs(finish = switch(family, gaussian = 2, student = 3, mixture = 4, poisson = 5))
    stopifnot(
      nchar(output$progress$html) > 0, nchar(output$step_description$html) > 0,
      nchar(output$conditional_description$html) > 0, nchar(output$step_table) > 0,
      nchar(output$model_equations$html) > 0, nchar(output$sample_summary) > 0,
      nchar(output$joint_plot$src) > 0, nchar(output$conditional_plot$src) > 0,
      nchar(output$diagnostic_plot$src) > 0, nchar(output$acf_plot$src) > 0
    )
  }
  session$setInputs(family = "mixture", separation = 2, weight = 0.35, rho = 0.2)
  stopifnot(frame() == 0, !playing())
  session$setInputs(play = 5)
  session$elapse(1001)
  stopifnot(frame() == 2, playing(), current()$coordinate == 2)
  session$setInputs(separation = 6)
  stopifnot(frame() == 0, !playing())
  session$setInputs(finish = 6)
  density_before <- target()$log_density(2, 2)
  session$setInputs(weight = 0.8)
  stopifnot(frame() == 0, !playing(), target()$log_density(2, 2) != density_before)
  session$setInputs(separation = 8, weight = 0.05, rho = -0.98, scale1 = 0.25, scale2 = 3)
  session$setInputs(finish = 7)
  stopifnot(nchar(output$joint_plot$src) > 0, nchar(output$conditional_plot$src) > 0,
            nchar(output$diagnostic_plot$src) > 0)
  session$setInputs(family = "poisson", total = 0, a = 1)
  stopifnot(frame() == 0, !playing())
  session$setInputs(finish = 8)
  stopifnot(nchar(output$joint_plot$src) > 0, nchar(output$conditional_plot$src) > 0)
  session$setInputs(n = 1, total = 200, hyper_rate = 0.1)
  session$setInputs(finish = 9)
  stopifnot(nchar(output$joint_plot$src) > 0, nchar(output$conditional_plot$src) > 0)
  session$setInputs(family = "gaussian", rho = -0.98, scale1 = 1, scale2 = 1)
  stopifnot(frame() == 0, !playing())
  session$setInputs(seed = 667)
  stopifnot(frame() == 0, !identical(chain()$states, original))
})
cat("Shiny plots, timed playback, pause/resume, stepping, replay, seek, resets, and diagnostics passed.\n")
