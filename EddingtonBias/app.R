source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Noise creates massive galaxies",
  "Watch symmetric errors reshape the rare, high-mass tail of an astronomical population.",
  "Equal chances of scattering up or down do not imply equal numbers moving in each direction.",
  controls = list(
    lesson_slider("noise", "Measurement error in log10 mass (dex)", 0, 0.6, 0.3, 0.02),
    lesson_slider("slope", "Population fall-off per dex", 0.3, 3, 1.5, 0.1),
    lesson_slider("threshold", "High-mass threshold: log10(M/Msun)", 9.5, 11.5, 10.5, 0.1)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/AstroStatistics_1.html",
  section = "Astrostatistics 1 - Eddington Bias and Monte Carlo Realisations",
  guide = c("Start with no measurement error, then increase it without changing the underlying galaxies.",
            "Steepen the population. Compare the number crossing the threshold in each direction.",
            "Replay: every galaxy receives a symmetric Gaussian error, but the high-mass counts need not stay unchanged."),
  explanation = c(
    "True log-masses follow a truncated exponential between 8.5 and 12 dex, representing a steeply declining toy mass function. Measurement errors are Gaussian and unbiased in log-mass, not necessarily in linear mass.",
    "Observed counts are a convolution of the true population with its error distribution. There are more lower-mass objects available to scatter upwards than higher-mass objects available to scatter downwards. The effect depends on the slope, curvature, noise, and boundaries; an individual realisation need not have more upward crossings.",
    "The teal curve is the analytic true density. The orange curve is the analytic convolution, including the finite population boundaries. Histograms use the same fixed bins and include all simulated measurements, not just selected ones."
  ),
  presets = list("No scatter" = list(noise = 0), "Steep and noisy" = list(noise = 0.5, slope = 2.5)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
