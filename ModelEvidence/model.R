line_evidence <- function(estimate, error, prior_sd) {
  log_bf <- dnorm(estimate, 0, sqrt(error^2 + prior_sd^2), log = TRUE) -
    dnorm(estimate, 0, error, log = TRUE)
  list(log_bf = log_bf, probability = plogis(log_bf),
       log_maximum_ratio = 0.5 * (estimate / error)^2)
}

calculate <- function(input, phase) {
  wavelength <- seq(500, 510, length.out = 70)
  template <- exp(-0.5 * ((wavelength - 505) / 0.35)^2)
  flux <- 1 + input$amplitude * template + input$noise * rnorm(length(template))
  estimate <- sum(template * (flux - 1)) / sum(template^2)
  error <- input$noise / sqrt(sum(template^2))
  list(wavelength = wavelength, template = template, flux = flux,
       estimate = estimate, error = error, evidence = line_evidence(estimate, error, input$prior))
}

draw <- function(result, input) {
  layout(matrix(1:2, ncol = 1), heights = c(1, 1.1))
  par(mar = c(3.8, 4.5, 2, 1), mgp = c(2.4, 0.8, 0))
  with(result, {
    plot(wavelength, flux, pch = 16, cex = 0.65, col = "#D5653499",
         xlab = "Wavelength (nm)", ylab = "Flux", main = "The same spectrum, two hypotheses")
    abline(h = 1, lty = 2, col = lesson_colours["truth"], lwd = 2)
    lines(wavelength, 1 + estimate * template, col = lesson_colours["inference"], lwd = 2)
    x <- seq(-max(3 * sqrt(error^2 + input$prior^2), abs(estimate) * 1.15),
             max(3 * sqrt(error^2 + input$prior^2), abs(estimate) * 1.15), length.out = 1400)
    null <- dnorm(x, 0, error)
    alternative <- dnorm(x, 0, sqrt(error^2 + input$prior^2))
    plot(x, null, type = "l", col = lesson_colours["truth"], lwd = 2,
         xlab = "Possible measured amplitude Ahat", ylab = "Predictive density",
         main = "Evidence scores predictions before the data")
    lines(x, alternative, col = lesson_colours["inference"], lwd = 2)
    abline(v = estimate, col = lesson_colours["data"], lwd = 2, lty = 2)
    legend("topright", c("No line: H0", "Proper amplitude prior: H1", "Measured amplitude"),
           col = c(lesson_colours["truth"], lesson_colours["inference"], lesson_colours["data"]),
           lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.75)
  })
}

measure <- function(result, input) {
  c("Measured amplitude / standard error" = sprintf("%.2f / %.2f", result$estimate, result$error),
    "log10 Bayes factor (line / no line)" = sprintf("%.2f", result$evidence$log_bf / log(10)),
    "P(line | data), equal prior model odds" = sprintf("%.1f%%", 100 * result$evidence$probability),
    "log10 maximum likelihood ratio" = sprintf("%.2f", result$evidence$log_maximum_ratio / log(10)))
}
