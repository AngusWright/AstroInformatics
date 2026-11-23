fit_colour_classifier <- function(x, label) {
  means <- rbind(colMeans(x[label == 0, , drop = FALSE]), colMeans(x[label == 1, , drop = FALSE]))
  pooled <- (cov(x[label == 0, ]) + cov(x[label == 1, ])) / 2
  weight <- solve(pooled, means[2, ] - means[1, ])
  list(weight = weight, threshold = sum(weight * colMeans(means)))
}

classify_colours <- function(model, x) {
  as.integer(drop(x %*% model$weight) > model$threshold)
}

calculate <- function(input, phase) {
  label <- rep(0:1, each = 200)
  make_sample <- function() {
    cbind((2 * label - 1) * 0.8 + 0.35 * rnorm(400),
          (2 * label - 1) * 0.4 + 0.45 * rnorm(400))
  }
  training <- make_sample()
  original <- make_sample()
  target <- make_sample() + matrix(c(input$offset, 0.4 * input$offset), 400, 2, byrow = TRUE) +
    input$noise * matrix(rnorm(800), 400, 2)
  model <- fit_colour_classifier(training, label)
  list(training = training, original = original, target = target, label = label, model = model,
       original_prediction = classify_colours(model, original),
       target_prediction = classify_colours(model, target))
}

draw <- function(result, input) {
  x <- if (input$view == "target") result$target else result$original
  prediction <- if (input$view == "target") result$target_prediction else result$original_prediction
  limits <- lapply(1:2, function(j) range(x[, j], result$training[, j]))
  gx <- seq(limits[[1]][1], limits[[1]][2], length.out = 120)
  gy <- seq(limits[[2]][1], limits[[2]][2], length.out = 120)
  background <- outer(gx, gy, function(a, b) {
    as.integer(a * result$model$weight[1] + b * result$model$weight[2] > result$model$threshold)
  })
  image(gx, gy, background, breaks = c(-0.5, 0.5, 1.5), col = c("#E9F4F3", "#EFEDF8"),
        xlab = "Measured colour 1", ylab = "Measured colour 2",
        main = if (input$view == "target") "A frozen classifier in a different survey" else "Independent validation under training conditions",
        useRaster = TRUE)
  points(x, pch = 21 + result$label, bg = lesson_colours[c("truth", "inference")][result$label + 1],
         col = "white", cex = 0.75)
  mistakes <- prediction != result$label
  points(x[mistakes, , drop = FALSE], pch = 1, cex = 1.15, col = lesson_colours["data"], lwd = 1.5)
  if (abs(result$model$weight[2]) > 1e-8) {
    abline(a = result$model$threshold / result$model$weight[2],
           b = -result$model$weight[1] / result$model$weight[2], lwd = 2, col = "#25344B")
  } else abline(v = result$model$threshold / result$model$weight[1], lwd = 2)
  legend("topleft", c("Star", "Quasar", "Classification error"),
         col = lesson_colours[c("truth", "inference", "data")], pch = c(16, 15, 1),
         bg = "#FFFFFFDD", cex = 0.8)
}

measure <- function(result, input) {
  c("Original-survey validation accuracy" = sprintf("%.1f%%", 100 * mean(result$original_prediction == result$label)),
    "New-survey validation accuracy" = sprintf("%.1f%%", 100 * mean(result$target_prediction == result$label)),
    "Training objects / validation objects per survey" = "400 / 400")
}
