calculate <- function(input, phase) {
  time <- 0
  arrivals <- numeric(0)
  repeat {
    time <- time + rexp(1, input$rate)
    if (time > 20) break
    arrivals <- c(arrivals, time)
  }
  exposure <- input$exposure * phase
  lambda <- input$rate * exposure
  counts <- qpois(runif(3000), lambda)
  list(arrivals = arrivals[arrivals <= exposure], exposure = exposure, lambda = lambda, counts = counts)
}

draw <- function(result, input) {
  layout(matrix(1:2, ncol = 1), heights = c(1, 1.6))
  par(mar = c(4.1, 4.5, 1.8, 1))
  plot(NA, xlim = c(0, input$exposure), ylim = c(0, 1), yaxt = "n", ylab = "",
       xlab = "Time (seconds)", main = "One possible photon stream")
  if (length(result$arrivals)) {
    segments(result$arrivals, 0.25, result$arrivals, 0.75, col = lesson_colours["data"], lwd = 2)
  } else {
    text(input$exposure / 2, 0.5, "No photons received yet", col = lesson_colours["text"])
  }
  abline(v = result$exposure, col = lesson_colours["inference"], lwd = 2, lty = 2)
  lower <- min(-3, floor(result$lambda - 4 * sqrt(result$lambda)))
  upper <- max(5, max(result$counts), ceiling(result$lambda + 4 * sqrt(result$lambda)))
  counts <- seq(0, upper)
  empirical <- tabulate(result$counts + 1, nbins = upper + 1) / length(result$counts)
  plot(counts, empirical, type = "h", lwd = 4, col = "#D5653455", xlim = c(lower, upper),
       ylim = c(0, max(empirical, dpois(counts, result$lambda)) * 1.1),
       xlab = "Photons counted in an exposure", ylab = "Probability per integer count",
       main = sprintf("Repeated-exposure counts: lambda = %.2f", result$lambda))
  points(counts, dpois(counts, result$lambda), pch = 16, cex = 0.6, col = lesson_colours["truth"])
  if (input$gaussian == "show" && result$lambda > 0) {
    x <- seq(lower, upper)
    mass <- pnorm(x + 0.5, result$lambda, sqrt(result$lambda)) -
      pnorm(x - 0.5, result$lambda, sqrt(result$lambda))
    lines(x, mass, col = lesson_colours["inference"], lwd = 2)
    bad <- x < 0
    segments(x[bad], 0, x[bad], mass[bad], col = "#D5653477", lwd = 4)
  }
  legend("topright", c("Simulated counts", "Exact Poisson PMF",
                       if (input$gaussian == "show" && result$lambda > 0) "Gaussian (continuity corrected)"),
         col = c(lesson_colours["data"], lesson_colours["truth"],
                 if (input$gaussian == "show" && result$lambda > 0) lesson_colours["inference"]),
         lty = c(1, NA, if (input$gaussian == "show" && result$lambda > 0) 1),
         pch = c(NA, 16, if (input$gaussian == "show" && result$lambda > 0) NA),
         bty = "n", cex = 0.7)
}

measure <- function(result, input) {
  c("Expected photons / photons in displayed stream" = sprintf("%.2f / %d", result$lambda, length(result$arrivals)),
    "Poisson count standard deviation" = sprintf("%.2f", sqrt(result$lambda)),
    "Relative count standard deviation" = if (result$lambda > 0) sprintf("%.1f%%", 100 / sqrt(result$lambda)) else "Undefined at zero exposure",
    "Gaussian mass below -0.5" = if (result$lambda > 0) sprintf("%.2f%%", 100 * pnorm(-0.5, result$lambda, sqrt(result$lambda))) else "No Gaussian at zero exposure")
}
