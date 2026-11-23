normal_hierarchy <- function(x, error, scatter, prior_mean = 0, prior_sd = 5) {
  variance <- 1 / (1 / prior_sd^2 + sum(1 / (scatter^2 + error^2)))
  mean <- variance * (prior_mean / prior_sd^2 + sum(x / (scatter^2 + error^2)))
  weight <- scatter^2 / (scatter^2 + error^2)
  list(mean = weight * x + (1 - weight) * mean,
       variance = scatter^2 * error^2 / (scatter^2 + error^2) + (1 - weight)^2 * variance,
       population_mean = mean, population_sd = sqrt(variance), weight = weight)
}

calculate <- function(input, phase) {
  truth <- 1.5 + input$scatter * rnorm(14)
  error <- seq(0.25, 1.7, length.out = 14) * input$noise
  observed <- truth + error * rnorm(14)
  posterior <- normal_hierarchy(observed, error, if (input$pooling == "complete") 0 else input$scatter)
  if (input$pooling == "none") {
    posterior$mean <- observed
    posterior$variance <- error^2
  }
  list(truth = truth, error = error, observed = observed, posterior = posterior,
       estimate = observed + phase * (posterior$mean - observed),
       sd = error + phase * (sqrt(posterior$variance) - error))
}

draw <- function(result, input) {
  with(result, {
    extent <- range(observed - 1.96 * error, observed + 1.96 * error, truth,
                    estimate - 1.96 * sd, estimate + 1.96 * sd)
    plot(NA, xlim = extent, ylim = c(0, 16), xlab = "Stellar velocity (km/s)",
         ylab = "", yaxt = "n", main = "Individual estimates, shared information")
    axis(2, at = 1:14, labels = paste("Star", 1:14), las = 1, cex.axis = 0.8)
    abline(v = 1.5, col = lesson_colours["truth"], lty = 2, lwd = 2)
    segments(observed - 1.96 * error, (1:14) + 0.12,
             observed + 1.96 * error, (1:14) + 0.12, col = "#B7C0CF")
    points(observed, (1:14) + 0.12, pch = 1, col = "#7F8A9B")
    segments(estimate - 1.96 * sd, (1:14) - 0.12,
             estimate + 1.96 * sd, (1:14) - 0.12, col = lesson_colours["inference"], lwd = 2)
    points(estimate, (1:14) - 0.12, pch = 16, col = lesson_colours["inference"])
    points(truth, 1:14, pch = 18, col = lesson_colours["truth"])
    legend("top", c("Measured +/- 1.96 SD", "Inferred +/- 1.96 SD", "Truth"),
           col = c("#7F8A9B", lesson_colours["inference"], lesson_colours["truth"]),
           pch = c(1, 16, 18), horiz = TRUE, bty = "n", cex = 0.75)
  })
}

measure <- function(result, input) {
  c("Median individual shift (km/s)" = sprintf("%.2f", median(abs(result$estimate - result$observed))),
    "True population mean (km/s)" = "1.50",
    "Assumed intrinsic dispersion" = if (input$pooling == "complete") "0 (all stars identical)" else if (input$pooling == "none") "Not pooled" else sprintf("%.1f km/s", input$scatter))
}
