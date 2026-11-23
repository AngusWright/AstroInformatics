robust_summaries <- function(x) {
  c(mean = mean(x), median = median(x), sd = sd(x),
    nmad = mad(x, constant = 1 / qnorm(0.75)))
}

calculate <- function(input, phase) {
  clean <- rnorm(120)[seq_len(input$n)]
  contaminants <- input$position + 0.1 * rnorm(10)
  observed <- clean
  indices <- seq_len(input$contaminants)
  observed[indices] <- contaminants[indices]
  list(clean = clean, observed = observed, indices = indices,
       jitter = runif(input$n, 2.2, 2.7), summaries = robust_summaries(observed))
}

draw <- function(result, input) {
  s <- result$summaries
  limits <- range(-3.5, 3.5, result$observed, s["mean"] + c(-1, 1) * s["sd"],
                  s["median"] + c(-1, 1) * s["nmad"])
  plot(NA, xlim = limits, ylim = c(-0.1, 3.6), yaxt = "n",
       xlab = "Velocity (km/s)", ylab = "", main = "A stable population, an unstable summary?")
  abline(v = 0, col = lesson_colours["truth"], lty = 2, lwd = 2)
  rect(-1, 2, 1, 2.9, col = "#087F8C15", border = NA)
  points(result$observed, result$jitter, pch = 16, col = "#7889A8", cex = 0.8)
  points(result$observed[result$indices], result$jitter[result$indices],
         pch = 16, col = lesson_colours["data"], cex = 1.1)
  segments(s["mean"] - s["sd"], 1.2, s["mean"] + s["sd"], 1.2,
           col = lesson_colours["data"], lwd = 4)
  points(s["mean"], 1.2, pch = 18, cex = 1.4, col = lesson_colours["data"])
  segments(s["median"] - s["nmad"], 0.45, s["median"] + s["nmad"], 0.45,
           col = lesson_colours["inference"], lwd = 4)
  points(s["median"], 0.45, pch = 18, cex = 1.4, col = lesson_colours["inference"])
  text(min(limits), c(3.2, 1.55, 0.8),
       c("Individual measured velocities", "Mean +/- sample SD", "Median +/- nMAD"),
       adj = 0, cex = 0.9)
  legend("topright", "True central population: mean 0, SD 1", lty = 2, lwd = 2,
         col = lesson_colours["truth"], bty = "n", cex = 0.8)
}

measure <- function(result, input) {
  setNames(sprintf("%.3f", result$summaries), c("Sample mean", "Sample median", "Sample SD", "Gaussian-scaled nMAD"))
}
