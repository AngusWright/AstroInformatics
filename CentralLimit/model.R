prepare <- function(input) {
  sample_population <- function(n) {
    switch(input$population,
           exponential = rexp(n),
           mixture = 5 * rbinom(n, 1, 0.05) + rnorm(n),
           cauchy = rcauchy(n))
  }
  list(values = matrix(sample_population(input$repetitions * 100), input$repetitions, 100),
       population = sample_population(8000))
}

calculate <- function(input, phase, prepared) {
  n <- max(1, round(1 + phase * (input$n - 1)))
  means <- rowMeans(prepared$values[, seq_len(n), drop = FALSE])
  finite <- input$population != "cauchy"
  mean <- if (input$population == "mixture") 0.25 else 1
  sd <- if (input$population == "mixture") sqrt(2.1875) else 1
  list(n = n, means = means, population = prepared$population, finite = finite,
       statistic = if (finite) sqrt(n) * (means - mean) / sd else means,
       mean_sd = if (finite) sd / sqrt(n) else NA_real_)
}

draw <- function(result, input) {
  layout(matrix(1:2, ncol = 1))
  par(mar = c(4.1, 4.5, 2, 1))
  limits <- switch(input$population, exponential = c(-0.2, 6), mixture = c(-4, 9), cauchy = c(-10, 10))
  histogram <- lesson_histogram(result$population, seq(limits[1], limits[2], length.out = 65))
  plot(histogram$mids, histogram$density, type = "h", lwd = 4, col = "#087F8C77",
       xlab = "Individual observation", ylab = "Density", main = "The original population does not become Gaussian")
  limits <- if (result$finite) c(-4, 6) else c(-10, 10)
  histogram <- lesson_histogram(result$statistic, seq(limits[1], limits[2], length.out = 60))
  x <- seq(limits[1], limits[2], length.out = 600)
  reference <- if (result$finite) dnorm(x) else dcauchy(x)
  plot(histogram$mids, histogram$density, type = "h", lwd = 4, col = "#D5653477",
       xlab = if (result$finite) "Standardised mean: sqrt(n) * (mean - mu) / sigma" else "Sample mean (no finite population mean or variance)",
       ylab = "Density", ylim = c(0, max(histogram$density, reference) * 1.05),
       main = sprintf("Independent sample means: n = %d", result$n))
  lines(x, reference, col = lesson_colours["inference"], lwd = 2)
  legend("topright", if (result$finite) "CLT reference: Normal(0,1)" else "Exact reference: Cauchy",
         lty = 1, lwd = 2, col = lesson_colours["inference"], bty = "n", cex = 0.8)
}

measure <- function(result, input) {
  c("Measurements per average at this frame" = format(result$n),
    "Independent averaging experiments" = format(length(result$means)),
    "Theoretical SD of the unstandardised mean" = if (result$finite) sprintf("%.3f", result$mean_sd) else "Does not exist",
    "Ordinary iid finite-variance CLT applies?" = if (result$finite) "Yes, asymptotically" else "No")
}
