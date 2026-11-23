calculate <- function(input, phase) {
  distance <- runif(1800)
  x <- -2 * distance + input$scatter * rnorm(1800)
  y <- -2 * distance + input$scatter * rnorm(1800)
  selected <- abs(distance - input$centre) <= input$width / 2
  list(distance = distance, x = x, y = y, selected = selected,
       full_correlation = cor(x, y),
       slice_correlation = if (sum(selected) >= 3) cor(x[selected], y[selected]) else NA_real_)
}

draw <- function(result, input) {
  palette <- colorRampPalette(c("#087F8C", "#B6DDD2", "#EB9A54"))(100)
  colours <- palette[pmin(100, floor(result$distance * 99) + 1)]
  plot(result$x, result$y, pch = 16, cex = 0.45, col = adjustcolor(colours, alpha.f = 0.25),
       xlab = "Log flux in band X (arbitrary zero point)",
       ylab = "Log flux in band Y (arbitrary zero point)", main = "A shared distance effect, not a direct physical relation")
  selected <- result$selected
  points(result$x[selected], result$y[selected], pch = 16, cex = 0.8, col = colours[selected])
  abline(lm(result$y ~ result$x), col = "#7F8A9B", lwd = 2, lty = 2)
  if (sum(selected) >= 3) abline(lm(result$y[selected] ~ result$x[selected]),
                                col = lesson_colours["inference"], lwd = 2)
  legend("topleft", c("Nearby objects", "Distant objects", "Whole-sample fit", "Highlighted slice fit"),
         col = c(palette[1], palette[100], "#7F8A9B", lesson_colours["inference"]),
         pch = c(16, 16, NA, NA), lty = c(NA, NA, 2, 1), bty = "n", cex = 0.8)
}

measure <- function(result, input) {
  c("Whole-sample Pearson correlation" = sprintf("%.2f", result$full_correlation),
    "Highlighted objects" = format(sum(result$selected)),
    "Distance-slice Pearson correlation" = if (is.finite(result$slice_correlation)) sprintf("%.2f", result$slice_correlation) else "Too few objects",
    "Direct X-to-Y relation in generating model" = "None")
}
