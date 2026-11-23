stellar_colour_density <- function(x) {
  0.6 * dnorm(x, -0.8, 0.35) + 0.4 * dnorm(x, 0.9, 0.35)
}

calculate <- function(input, phase) {
  component <- runif(1000) < 0.6
  values <- ifelse(component, -0.8, 0.9) + 0.35 * rnorm(1000)
  values <- values[seq_len(input$n)]
  offset <- input$origin * input$width
  lower <- floor((min(values) - offset) / input$width) * input$width + offset
  upper <- ceiling((max(values) - offset) / input$width) * input$width + offset + input$width
  breaks <- seq(lower, upper, by = input$width)
  kde <- density(values, bw = input$bandwidth, from = min(-3, values), to = max(3, values), n = 512)
  peaks <- which(diff(sign(diff(kde$y))) == -2) + 1
  peaks <- peaks[kde$y[peaks] > 0.02 * max(kde$y)]
  list(values = values, histogram = hist(values, breaks = breaks, plot = FALSE),
       kde = kde, peaks = peaks)
}

draw <- function(result, input) {
  histogram <- result$histogram
  x <- result$kde$x
  plot(histogram, freq = FALSE, col = "#CBD2DF88", border = "white",
       xlim = range(x), ylim = c(0, max(histogram$density, result$kde$y, stellar_colour_density(x)) * 1.1),
       xlab = "Stellar colour (illustrative units)", ylab = "Probability density",
       main = "Same objects, different apparent structure")
  lines(x, stellar_colour_density(x), col = lesson_colours["truth"], lwd = 3, lty = 2)
  lines(x, result$kde$y, col = lesson_colours["inference"], lwd = 3)
  rug(result$values, col = "#D5653444")
  legend("topright", c("True population", "Estimated KDE", "Histogram"),
         col = c(lesson_colours["truth"], lesson_colours["inference"], "#8795AE"),
         lty = c(2, 1, 1), lwd = c(2, 2, 5), bty = "n", cex = 0.85)
}

measure <- function(result, input) {
  c("Objects (unchanged by smoothing controls)" = format(length(result$values)),
    "Apparent KDE peaks above display threshold" = format(length(result$peaks)),
    "Known generating mixture components" = "2",
    "Histogram area" = sprintf("%.3f", sum(result$histogram$density * diff(result$histogram$breaks))))
}
