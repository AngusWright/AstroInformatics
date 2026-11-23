selection_observed <- function(log_luminosity, distance, limit) {
  apparent <- -18 - 2.5 * log_luminosity + 5 * log10(distance) + 25
  apparent <= limit
}

calculate <- function(input, phase) {
  distance <- (30^3 + runif(2400) * (800^3 - 30^3))^(1 / 3)
  luminosity <- rnorm(2400, 0, 0.6)
  limit <- 26 + phase * (input$limit - 26)
  detected <- selection_observed(luminosity, distance, limit)
  list(distance = distance, luminosity = luminosity, limit = limit, detected = detected)
}

draw <- function(result, input) {
  with(result, {
    plot(distance, luminosity, type = "n", xlim = c(30, 800), ylim = c(-2.2, 2.4),
         xlab = "Distance (Mpc)", ylab = "Intrinsic log10 luminosity (L/L0)",
         main = sprintf("What survives a magnitude limit of %.1f?", limit))
    if (input$overlay == "show") points(distance[!detected], luminosity[!detected],
                                       pch = 16, cex = 0.5, col = "#B7C0CF70")
    points(distance[detected], luminosity[detected], pch = 16, cex = 0.6, col = "#D5653499")
    x <- seq(30, 800, length.out = 300)
    lines(x, (5 * log10(x) + 7 - limit) / 2.5, col = "#25344B", lwd = 2)
    abline(h = 0, col = lesson_colours["truth"], lty = 2, lwd = 2)
    bins <- cut(distance, seq(0, 800, by = 80), labels = FALSE)
    means <- vapply(1:10, function(i) {
      indices <- bins == i & detected
      if (any(indices)) mean(luminosity[indices]) else NA_real_
    }, numeric(1))
    lines(seq(40, 760, by = 80), means, col = lesson_colours["data"], lwd = 3)
    legend("topleft", c("Detected galaxies / binned mean", "True population mean", "Detection boundary"),
           col = c(lesson_colours["data"], lesson_colours["truth"], "#25344B"),
           lty = c(1, 2, 1), lwd = 2, bty = "n", cex = 0.85)
  })
}

measure <- function(result, input) {
  c("Galaxies detected" = sprintf("%d / %d", sum(result$detected), length(result$detected)),
    "Completeness in this realisation" = sprintf("%.1f%%", 100 * mean(result$detected)),
    "Detected mean log-luminosity" = if (any(result$detected)) {
      sprintf("%.2f (truth: 0)", mean(result$luminosity[result$detected]))
    } else "No detections")
}
