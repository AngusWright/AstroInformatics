source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Knowing the mean is not predicting the next star",
  "More data can pin down a population mean while leaving substantial uncertainty about its next member.",
  "Posterior parameter uncertainty and posterior predictive uncertainty are different quantities.",
  controls = list(
    lesson_slider("n", "Number of measured member stars", 1, 300, 20, 1),
    lesson_slider("intrinsic", "Intrinsic velocity dispersion (km/s)", 0, 4, 1.5, 0.1),
    lesson_slider("measurement", "Measurement standard deviation (km/s)", 0.1, 3, 0.8, 0.1)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section3c_html.html",
  section = "3c - Predictive Distributions",
  guide = c("Increase the sample size. Which interval keeps shrinking?",
            "Set intrinsic dispersion to zero. A future measurement still has measurement noise.",
            "Increase intrinsic dispersion while holding the measurement precision fixed."),
  explanation = c(
    "Member velocities follow Normal(mu,tau^2); measurements add independent Normal(0,sigma^2) noise. Both tau and sigma are known. The mean has proper prior mu ~ Normal(0,5^2); its simulated true value is 2 km/s.",
    "The mean posterior variance is V = 1/[1/25+n/(tau^2+sigma^2)]. Its mean is V*sum(x)/(tau^2+sigma^2). The latent velocity of a new star has predictive variance V+tau^2; a new measured velocity has variance V+tau^2+sigma^2.",
    "Each row is a probability density scaled to its own peak height, to compare widths without hiding a broad prediction under a narrow posterior. Horizontal bars are central 95% intervals. They are credible or predictive intervals according to the row, not frequentist confidence intervals.",
    "Increasing n reuses a nested set of simulated stars. Selection, membership uncertainty, and uncertainty in the dispersions are excluded; including those would require a richer model."
  ),
  presets = list("Few stars" = list(n = 2), "Many stars" = list(n = 300),
                 "No intrinsic scatter" = list(intrinsic = 0))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
