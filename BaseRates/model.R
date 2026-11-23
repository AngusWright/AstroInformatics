positive_predictive_value <- function(prevalence, sensitivity, false_positive) {
  prevalence * sensitivity / (prevalence * sensitivity + (1 - prevalence) * false_positive)
}

calculate <- function(input, phase) {
  prevalence <- 10^input$log_prevalence
  truth <- runif(10000) < prevalence
  flagged <- runif(10000) < ifelse(truth, input$sensitivity, input$false_positive)
  scanned <- seq_along(truth) <= floor(10000 * phase)
  indices <- which(flagged & scanned)
  list(prevalence = prevalence, truth = truth, indices = indices,
       genuine = truth[indices], scanned = sum(scanned),
       probability = positive_predictive_value(prevalence, input$sensitivity, input$false_positive))
}

draw <- function(result, input) {
  count <- length(result$indices)
  columns <- max(5, ceiling(sqrt(max(1, count) * 1.6)))
  rows <- max(4, ceiling(max(1, count) / columns))
  plot(NA, xlim = c(0.5, columns + 0.5), ylim = c(rows + 1.5, 0.5),
       xlab = "", ylab = "", axes = FALSE, asp = 1,
       main = sprintf("Flagged candidates after scanning %s objects", format(result$scanned, big.mark = ",")))
  if (count) {
    colours <- if (input$reveal == "show") {
      ifelse(result$genuine, lesson_colours["truth"], lesson_colours["data"])
    } else rep(lesson_colours["inference"], count)
    points((seq_len(count) - 1) %% columns + 1,
           (seq_len(count) - 1) %/% columns + 1, pch = 15, col = colours, cex = 0.9)
  } else {
    text((columns + 1) / 2, rows / 2, "No candidates flagged yet", col = lesson_colours["text"])
  }
  legend("bottom", if (input$reveal == "show") c("Real event", "False candidate") else "Flagged; identity unknown",
         col = if (input$reveal == "show") lesson_colours[c("truth", "data")] else lesson_colours["inference"],
         pch = 15, horiz = TRUE, bty = "n", cex = 0.9)
}

measure <- function(result, input) {
  c("Population fraction of real events" = sprintf("%.3g%%", 100 * result$prevalence),
    "Detection efficiency (sensitivity)" = sprintf("%.0f%%", 100 * input$sensitivity),
    "Theoretical P(real | flag)" = sprintf("%.1f%%", 100 * result$probability),
    "Genuine candidates / all flags in this run" = sprintf("%d / %d", sum(result$genuine), length(result$indices)))
}
