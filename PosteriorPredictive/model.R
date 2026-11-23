predictive_normal <- function(x, intrinsic, measurement, prior_sd = 5) {
  observation_variance <- intrinsic^2 + measurement^2
  variance <- 1 / (1 / prior_sd^2 + length(x) / observation_variance)
  list(mean = variance * sum(x) / observation_variance,
       variance = c(mean = variance, star = variance + intrinsic^2,
                    measurement = variance + intrinsic^2 + measurement^2))
}

calculate <- function(input, phase) {
  stars <- 2 + input$intrinsic * rnorm(300)
  observations <- stars + input$measurement * rnorm(300)
  x <- observations[seq_len(input$n)]
  list(x = x, posterior = predictive_normal(x, input$intrinsic, input$measurement))
}

draw <- function(result, input) {
  centre <- result$posterior$mean
  sd <- sqrt(result$posterior$variance)
  x <- seq(centre - 4 * max(sd), centre + 4 * max(sd), length.out = 1000)
  labels <- c("Population mean", "Next true stellar velocity", "Next measured velocity")
  colours <- lesson_colours[c("inference", "truth", "data")]
  plot(NA, xlim = range(x, 2), ylim = c(0.5, 4.2), yaxt = "n",
       xlab = "Velocity (km/s)", ylab = "", main = "Three different questions, three different uncertainties")
  abline(v = 2, col = "#85909F", lty = 3)
  for (j in 1:3) {
    base <- 4 - j
    grid <- sort(unique(c(x, seq(centre - 4 * sd[j], centre + 4 * sd[j], length.out = 400))))
    relative <- dnorm(grid, centre, sd[j]) / dnorm(centre, centre, sd[j]) * 0.65
    polygon(c(grid, rev(grid)), c(base + relative, rep(base, length(grid))),
            col = adjustcolor(colours[j], alpha.f = 0.2), border = NA)
    lines(grid, base + relative, col = colours[j], lwd = 2)
    segments(centre - 1.96 * sd[j], base - 0.08, centre + 1.96 * sd[j], base - 0.08,
             lwd = 4, col = colours[j])
    text(min(x), base + 0.78, labels[j], adj = 0, cex = 0.85, col = colours[j])
  }
  mtext("Density rows are independently peak-scaled; bars enclose 95% probability.", side = 3, line = 0.2, cex = 0.75)
}

measure <- function(result, input) {
  sd <- sqrt(result$posterior$variance)
  c("Posterior SD of population mean" = sprintf("%.3f km/s", sd[1]),
    "Predictive SD of next true star" = sprintf("%.3f km/s", sd[2]),
    "Predictive SD of next measurement" = sprintf("%.3f km/s", sd[3]))
}
