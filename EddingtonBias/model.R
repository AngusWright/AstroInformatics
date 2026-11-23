mass_density <- function(x, slope, noise = 0, lower = 8.5, upper = 12) {
  normalisation <- 1 - exp(-slope * (upper - lower))
  if (noise == 0) {
    return(ifelse(x >= lower & x <= upper, slope * exp(-slope * (x - lower)) / normalisation, 0))
  }
  slope / normalisation * exp(slope * (lower - x) + (slope * noise)^2 / 2) *
    (pnorm((upper - x + slope * noise^2) / noise) -
       pnorm((lower - x + slope * noise^2) / noise))
}

calculate <- function(input, phase) {
  true <- 8.5 - log1p(-runif(8000) * (1 - exp(-input$slope * 3.5))) / input$slope
  noise <- phase * input$noise
  observed <- true + noise * rnorm(length(true))
  list(true = true, observed = observed, noise = noise,
       up = sum(true < input$threshold & observed >= input$threshold),
       down = sum(true >= input$threshold & observed < input$threshold))
}

draw <- function(result, input) {
  breaks <- seq(floor(min(result$true, result$observed) * 10) / 10,
                ceiling(max(result$true, result$observed) * 10) / 10 + 0.1, by = 0.1)
  truth <- hist(result$true, breaks = breaks, plot = FALSE)
  observed <- hist(result$observed, breaks = breaks, plot = FALSE)
  x <- seq(7, 13.5, length.out = 600)
  plot(observed$mids, observed$density, type = "n", xlim = c(7.5, 13),
       ylim = c(0, max(truth$density, observed$density) * 1.1),
       xlab = "log10 stellar mass (M/Msun)", ylab = "Probability density",
       main = "Unbiased errors; biased population counts")
  polygon(c(observed$mids, rev(observed$mids)),
          c(observed$density, rep(0, length(observed$mids))),
          col = "#D5653430", border = NA)
  lines(x, mass_density(x, input$slope), col = lesson_colours["truth"], lwd = 3)
  lines(x, mass_density(x, input$slope, result$noise), col = lesson_colours["data"], lwd = 3)
  abline(v = input$threshold, lty = 2, lwd = 2)
  legend("topright", c("True population", "Observed population", "High-mass threshold"),
         col = c(lesson_colours["truth"], lesson_colours["data"], lesson_colours["text"]),
         lty = c(1, 1, 2), lwd = 2, bty = "n")
}

measure <- function(result, input) {
  c("Scattered into high-mass sample" = format(result$up),
    "Scattered out of high-mass sample" = format(result$down),
    "Net change in high-mass counts" = sprintf("%+d", result$up - result$down))
}
