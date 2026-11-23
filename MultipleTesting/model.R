search_false_probability <- function(threshold, trials, two_sided = FALSE) {
  probability <- pnorm(threshold, lower.tail = FALSE) * if (two_sided) 2 else 1
  list(single = probability, global = -expm1(trials * log1p(-probability)),
       expected = trials * probability)
}

calculate <- function(input, phase) {
  side <- input$side
  z <- rnorm(side^2)
  tested <- floor(phase * length(z))
  tested_mask <- seq_along(z) <= tested
  hit <- (if (input$tails == "two") abs(z) else z) >= input$threshold & tested_mask
  list(side = side, z = z, tested = tested, tested_mask = tested_mask, hit = hit,
       probability = search_false_probability(input$threshold, tested, input$tails == "two"))
}

draw <- function(result, input) {
  with(result, {
    image_values <- matrix(ifelse(tested_mask, pmax(-3, pmin(3, z)), NA_real_), side, side)
    palette <- colorRampPalette(c("#263956", "#E9EEF4", "#F8B991"))(100)
    image(seq_len(side), seq_len(side), image_values, col = palette, zlim = c(-3, 3),
          xlab = "Independent sky resolution element", ylab = "", asp = 1,
          main = sprintf("Noise-only sky: %d locations searched", tested), useRaster = TRUE)
    indices <- which(hit)
    if (length(indices)) points((indices - 1) %% side + 1, (indices - 1) %/% side + 1,
                                pch = 1, cex = 1.3, col = "#D54524", lwd = 2)
    legend("topleft", "Circled peaks are all false positives", pch = 1, col = "#D54524",
           bg = "#FFFFFFDD", cex = 0.8, bty = "o")
  })
}

measure <- function(result, input) {
  c("False detections in this realisation" = format(sum(result$hit)),
    "Single-location false-positive probability" = sprintf("%.4g", result$probability$single),
    "Expected false detections in searched area" = sprintf("%.2f", result$probability$expected),
    "Chance of at least one in searched area" = sprintf("%.1f%%", 100 * result$probability$global))
}
