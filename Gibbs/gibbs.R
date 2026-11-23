gibbs_normal <- function(location, scale) {
  list(
    label = sprintf("Normal(mean = %.4g, SD = %.4g)", location, scale),
    draw = function() rnorm(1, location, scale),
    density = function(x) dnorm(x, location, scale),
    quantile = function(p) qnorm(p, location, scale)
  )
}

gibbs_student <- function(location, scale, df) {
  list(
    label = sprintf("Student-t(df = %.4g, location = %.4g, scale = %.4g)", df, location, scale),
    draw = function() location + scale * rt(1, df),
    density = function(x) dt((x - location) / scale, df) / scale,
    quantile = function(p) location + scale * qt(p, df)
  )
}

gibbs_gamma <- function(shape, rate) {
  list(
    label = sprintf("Gamma(shape = %.4g, rate = %.4g)", shape, rate),
    draw = function() rgamma(1, shape = shape, rate = rate),
    density = function(x) dgamma(x, shape = shape, rate = rate),
    quantile = function(p) qgamma(p, shape = shape, rate = rate)
  )
}

gibbs_target <- function(family = "gaussian", location = c(0, 0), scale = c(1, 1),
                         rho = 0.8, df = 5, n = 10, total = 30,
                         a = 2, c = 2, d = 1) {
  stopifnot(length(family) == 1, family %in% c("gaussian", "student", "poisson"))
  if (family %in% c("gaussian", "student")) {
    stopifnot(length(location) == 2, all(is.finite(location)),
              length(scale) == 2, all(is.finite(scale)), all(scale > 0),
              length(rho) == 1, is.finite(rho), abs(rho) < 1)
    if (family == "student") stopifnot(length(df) == 1, is.finite(df), df > 2)
    conditional <- function(state, coordinate) {
      other <- 3L - coordinate
      residual <- (state[other] - location[other]) / scale[other]
      centre <- location[coordinate] + rho * scale[coordinate] * residual
      width <- scale[coordinate] * sqrt(1 - rho^2)
      if (family == "student") {
        return(gibbs_student(centre, width * sqrt((df + residual^2) / (df + 1)), df + 1))
      }
      gibbs_normal(centre, width)
    }
    log_density <- function(x, y) {
      u <- (x - location[1]) / scale[1]
      v <- (y - location[2]) / scale[2]
      q <- (u^2 - 2 * rho * u * v + v^2) / (1 - rho^2)
      if (family == "student") {
        lgamma((df + 2) / 2) - lgamma(df / 2) - log(df * pi) -
          sum(log(scale)) - log1p(-rho^2) / 2 - (df + 2) / 2 * log1p(q / df)
      } else {
        -log(2 * pi) - sum(log(scale)) - log1p(-rho^2) / 2 - q / 2
      }
    }
    marginal <- function(x, coordinate) {
      if (family == "student") {
        dt((x - location[coordinate]) / scale[coordinate], df) / scale[coordinate]
      } else {
        dnorm(x, location[coordinate], scale[coordinate])
      }
    }
    quantiles <- if (family == "student") qt(c(0.002, 0.998), df) else qnorm(c(0.001, 0.999))
    limits <- lapply(1:2, function(j) location[j] + scale[j] * quantiles)
    return(list(
      family = family,
      label = if (family == "student") "Bivariate Student-t posterior" else "Bivariate Gaussian posterior",
      axes = c("theta1", "theta2"), positive = FALSE,
      conditional = conditional, log_density = log_density, marginal = marginal,
      limits = limits, mean = location,
      sd = scale * if (family == "student") sqrt(df / (df - 2)) else 1
    ))
  }
  stopifnot(
    length(n) == 1, is.finite(n), n >= 1, n == round(n),
    length(total) == 1, is.finite(total), total >= 0, total == round(total),
    length(a) == 1, is.finite(a), a >= 1,
    length(c) == 1, is.finite(c), c >= 1,
    length(d) == 1, is.finite(d), d > 0
  )
  list(
    family = family, label = "Poisson likelihood with a Gamma-Gamma hierarchy",
    axes = c("lambda (theta1)", "beta (theta2)"), positive = TRUE,
    conditional = function(state, coordinate) {
      if (coordinate == 1L) gibbs_gamma(a + total, n + state[2]) else gibbs_gamma(a + c, d + state[1])
    },
    log_density = function(x, y) {
      (total + a - 1) * log(x) + (a + c - 1) * log(y) - n * x - y * x - d * y
    },
    limits = NULL,
    marginal = NULL, mean = NULL, sd = NULL
  )
}

gibbs_chain <- function(target, start, sweeps = 100, seed = 666) {
  stopifnot(
    length(start) == 2, all(is.finite(start)),
    length(sweeps) == 1, is.finite(sweeps), sweeps >= 1, sweeps == round(sweeps),
    length(seed) == 1, is.finite(seed), seed >= 0, seed <= .Machine$integer.max,
    seed == round(seed)
  )
  if (target$positive) stopifnot(all(start > 0))
  previous_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  on.exit({
    if (is.null(previous_seed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  states <- matrix(NA_real_, nrow = 2 * sweeps + 1, ncol = 2,
                   dimnames = list(NULL, target$axes))
  states[1, ] <- start
  for (update in seq_len(2 * sweeps)) {
    coordinate <- (update - 1L) %% 2L + 1L
    states[update + 1L, ] <- states[update, ]
    distribution <- target$conditional(states[update, ], coordinate)
    states[update + 1L, coordinate] <- distribution$draw()
  }
  list(target = target, states = states, sweeps = sweeps, seed = seed)
}

gibbs_frame <- function(chain, update) {
  stopifnot(length(update) == 1, is.finite(update), update >= 0,
            update <= 2 * chain$sweeps, update == round(update))
  coordinate <- if (update == 0) 1L else (update - 1L) %% 2L + 1L
  before <- chain$states[if (update == 0) 1L else update, ]
  after <- chain$states[update + 1L, ]
  list(
    update = update, coordinate = coordinate, before = before, after = after,
    distribution = chain$target$conditional(before, coordinate),
    has_draw = update > 0
  )
}

gibbs_samples <- function(chain, update, burn = 0) {
  stopifnot(length(burn) == 1, is.finite(burn), burn >= 0, burn == round(burn))
  gibbs_frame(chain, update)
  completed <- floor(update / 2)
  indices <- if (completed > burn) 2L * seq.int(burn + 1L, completed) + 1L else integer(0)
  chain$states[indices, , drop = FALSE]
}
