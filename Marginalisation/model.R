mixture_joint <- function(x, y, weight, width) {
  component <- function(x, y, mx, my, sx, sy) {
    u <- (x - mx) / sx
    v <- (y - my) / sy
    exp(-(u^2 - 0.9 * u * v + v^2) / (2 * (1 - 0.45^2))) /
      (2 * pi * sx * sy * sqrt(1 - 0.45^2))
  }
  (1 - weight) * component(x, y, -1, 0, 0.23, 0.3) +
    weight * component(x, y, 1.3, 1, width, width)
}

mixture_marginal <- function(x, weight, width) {
  (1 - weight) * dnorm(x, -1, 0.23) + weight * dnorm(x, 1.3, width)
}

mixture_interval <- function(lower, upper, weight, width) {
  (1 - weight) * (pnorm(upper, -1, 0.23) - pnorm(lower, -1, 0.23)) +
    weight * (pnorm(upper, 1.3, width) - pnorm(lower, 1.3, width))
}

calculate <- function(input, phase) {
  x <- seq(-5, 7, length.out = 250)
  y <- seq(-6, 8, length.out = 250)
  list(x = x, y = y, joint = outer(x, y, mixture_joint, input$weight, input$width),
       marginal = mixture_marginal(x, input$weight, input$width),
       lower = input$centre - 0.3, upper = input$centre + 0.3,
       probability = mixture_interval(input$centre - 0.3, input$centre + 0.3,
                                      input$weight, input$width))
}

draw <- function(result, input) {
  layout(matrix(1:2, ncol = 1), heights = c(1.5, 1))
  par(mar = c(4.1, 4.5, 2, 1))
  with(result, {
    image(x, y, sqrt(joint), col = colorRampPalette(c("#FFFFFF", "#CED5E9", "#5059A6"))(120),
          xlab = "Parameter x", ylab = "Nuisance parameter y",
          main = "Joint density: tall peaks are not probability masses", useRaster = TRUE)
    contour(x, y, joint, add = TRUE, drawlabels = FALSE, col = "#66749977")
    rect(lower, min(y), upper, max(y), col = "#087F8C30", border = lesson_colours["truth"], lwd = 2)
    plot(x, marginal, type = "l", lwd = 2, col = lesson_colours["inference"],
         xlab = "Parameter x", ylab = "Marginal density", main = "Integrate over all y")
    band <- seq(lower, upper, length.out = 150)
    polygon(c(lower, band, upper), c(0, mixture_marginal(band, input$weight, input$width), 0),
            col = "#087F8C66", border = NA)
    lines(x, marginal, col = lesson_colours["inference"], lwd = 2)
  })
}

measure <- function(result, input) {
  c("Probability in broad component" = sprintf("%.0f%%", 100 * input$weight),
    "Probability in highlighted x interval" = sprintf("%.1f%%", 100 * result$probability),
    "Approximate marginal mode (not a probability)" = sprintf("%.2f", result$x[which.max(result$marginal)]))
}
