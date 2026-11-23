positive_distribution <- function(family, slope) {
  if (family == "exponential") return(list(cdf = pexp, quantile = qexp, density = dexp))
  if (family == "lognormal") return(list(
    cdf = function(x) plnorm(x, 0, 0.7),
    quantile = function(p) qlnorm(p, 0, 0.7),
    density = function(x) dlnorm(x, 0, 0.7)))
  lower <- 0.1
  upper <- 20
  if (abs(slope - 1) < 1e-8) return(list(
    cdf = function(x) pmin(1, pmax(0, log(x / lower) / log(upper / lower))),
    quantile = function(p) lower * (upper / lower)^p,
    density = function(x) ifelse(x >= lower & x <= upper, 1 / (x * log(upper / lower)), 0)))
  exponent <- 1 - slope
  normalisation <- (upper^exponent - lower^exponent) / exponent
  list(
    cdf = function(x) pmin(1, pmax(0, (x^exponent - lower^exponent) / (upper^exponent - lower^exponent))),
    quantile = function(p) (lower^exponent + p * (upper^exponent - lower^exponent))^(1 / exponent),
    density = function(x) ifelse(x >= lower & x <= upper, x^(-slope) / normalisation, 0)
  )
}

calculate <- function(input, phase) {
  distribution <- positive_distribution(input$family, input$slope)
  uniforms <- runif(input$n)
  all <- distribution$quantile(uniforms)
  count <- floor(phase * input$n)
  revealed <- seq_len(count)
  probabilities <- sort(uniforms[revealed])
  discrepancy <- if (count) max(seq_len(count) / count - probabilities,
                                probabilities - (seq_len(count) - 1) / count) else NA_real_
  list(distribution = distribution, u = uniforms[revealed], values = all[revealed],
       limits = range(log10(distribution$quantile(c(0.00001, 0.99999))), log10(all)),
       discrepancy = discrepancy)
}

draw <- function(result, input) {
  layout(matrix(1:2, ncol = 1))
  par(mar = c(4.1, 4.5, 2, 1))
  z <- seq(result$limits[1], result$limits[2], length.out = 700)
  plot(z, result$distribution$cdf(10^z), type = "l", lwd = 3,
       col = lesson_colours["inference"], ylim = c(0, 1),
       xlab = "log10 target value", ylab = "Cumulative probability",
       main = "Uniform probability -> inverse CDF -> target value")
  indices <- tail(seq_along(result$u), 12)
  for (j in indices) {
    value <- log10(result$values[j])
    segments(min(z), result$u[j], value, result$u[j], col = "#D5653488")
    segments(value, result$u[j], value, 0, col = "#D5653488")
    points(value, result$u[j], pch = 16, cex = 0.6, col = lesson_colours["data"])
  }
  density <- log(10) * 10^z * result$distribution$density(10^z)
  histogram <- lesson_histogram(log10(result$values), seq(min(z), max(z), length.out = 45))
  plot(histogram$mids, histogram$density, type = "h", lwd = 5, col = "#D5653466",
       xlab = "log10 sampled value", ylab = "Density per log10 unit",
       ylim = c(0, max(histogram$density, density) * 1.1), main = "The transformed population")
  lines(z, density, col = lesson_colours["truth"], lwd = 2)
}

measure <- function(result, input) {
  c("Uniform draws transformed so far" = format(length(result$values)),
    "Maximum empirical-CDF discrepancy" = if (is.finite(result$discrepancy)) sprintf("%.3f", result$discrepancy) else "No draws yet",
    "Display scale" = "log10 value; density includes Jacobian")
}
