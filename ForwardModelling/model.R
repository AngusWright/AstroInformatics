restricted_encounter <- function(radius, angle, impact, speed, steps = 80, dt = 0.1) {
  position <- cbind(radius * cos(angle), radius * sin(angle))
  circular <- radius / (radius^2 + 0.3^2)^(3 / 4)
  velocity <- cbind(-circular * sin(angle), circular * cos(angle))
  acceleration <- function(position, time) {
    perturber <- c(-5 + speed * time, impact)
    separation <- matrix(perturber, nrow(position), 2, byrow = TRUE) - position
    -position / (rowSums(position^2) + 0.3^2)^(3 / 2) +
      0.5 * separation / (rowSums(separation^2) + 0.3^2)^(3 / 2) -
      matrix(0.5 * perturber / (sum(perturber^2) + 0.3^2)^(3 / 2),
             nrow(position), 2, byrow = TRUE)
  }
  trajectory <- array(NA_real_, c(steps + 1, length(radius), 2))
  trajectory[1, , ] <- position
  for (i in seq_len(steps)) {
    velocity <- velocity + dt / 2 * acceleration(position, (i - 1) * dt)
    position <- position + dt * velocity
    velocity <- velocity + dt / 2 * acceleration(position, i * dt)
    trajectory[i + 1, , ] <- position
  }
  trajectory
}

prepare <- function(input) {
  radius <- seq(0.65, 1.5, length.out = 40)
  angle <- runif(40, 0, 2 * pi)
  log_sd <- sqrt(log1p(input$uncertainty^2))
  impacts <- input$impact * exp(log_sd * rnorm(12) - log_sd^2 / 2)
  speeds <- input$speed * exp(log_sd * rnorm(12) - log_sd^2 / 2)
  list(nominal = restricted_encounter(radius, angle, input$impact, input$speed),
       ensemble = lapply(1:12, function(j) restricted_encounter(radius, angle, impacts[j], speeds[j])))
}

calculate <- function(input, phase, prepared) {
  index <- round(80 * phase) + 1
  nominal <- prepared$nominal[index, , ]
  ensemble <- lapply(prepared$ensemble, function(path) path[index, , ])
  outside <- vapply(ensemble, function(position) mean(rowSums(position^2) > 2.5^2), numeric(1))
  list(nominal = nominal, ensemble = ensemble, outside = outside, time = (index - 1) * 0.1,
       perturber = c(-5 + input$speed * (index - 1) * 0.1, input$impact))
}

draw <- function(result, input) {
  plot(NA, xlim = c(-5, 7), ylim = c(-5, 5), asp = 1,
       xlab = "Host-centred x (model units)", ylab = "Host-centred y (model units)",
       main = sprintf("A family of predictions at time %.1f", result$time))
  for (position in result$ensemble) points(position, pch = 16, cex = 0.55, col = "#5059A628")
  points(result$nominal, pch = 16, cex = 0.8, col = lesson_colours["truth"])
  points(0, 0, pch = 8, cex = 1.5, lwd = 2, col = lesson_colours["text"])
  points(result$perturber[1], result$perturber[2], pch = 8, cex = 1.5, lwd = 2, col = lesson_colours["data"])
  symbols(0, 0, circles = 2.5, inches = FALSE, add = TRUE, fg = "#B7C0CF")
  legend("topleft", c("Nominal encounter", "Uncertain encounter ensemble", "Nominal perturber"),
         col = lesson_colours[c("truth", "inference", "data")], pch = c(16, 16, 8),
         bty = "n", cex = 0.8)
  mtext("Finite viewport; particles outside it still enter the numerical summaries.", side = 3, line = 0.2, cex = 0.7)
}

measure <- function(result, input) {
  c("Model time" = sprintf("%.1f", result$time),
    "Ensemble members" = "12",
    "Beyond radius 2.5: median ensemble fraction" = sprintf("%.0f%%", 100 * median(result$outside)),
    "Beyond radius 2.5: ensemble range" = sprintf("%.0f to %.0f%%", 100 * min(result$outside), 100 * max(result$outside)))
}
