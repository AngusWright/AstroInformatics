source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Is that spectral line real?",
  "A flexible model can fit well while making a poor prediction before seeing the data.",
  "Bayesian evidence averages likelihood over a proper prior; it is not the maximum likelihood.",
  controls = list(
    lesson_slider("amplitude", "True line amplitude (positive or negative)", -2, 3, 0.5, 0.1),
    lesson_slider("noise", "Known flux uncertainty", 0.3, 2, 1, 0.1),
    lesson_slider("prior", "Alternative amplitude prior SD", 0.1, 10, 2, 0.1)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section4b_html.html",
  section = "4b - Bayesian Model Comparison and Additional Parameters",
  guide = c("Hold the spectrum fixed and widen the alternative prior. Does the best fit change? Does its evidence?",
            "Compare the noise-only and clear-line examples.",
            "In the lower panel, compare the two predictive densities exactly at the measured amplitude."),
  explanation = c(
    "The continuum, line centre, and line width are known. Seventy independent Gaussian flux measurements constrain only the line amplitude A. Positive A is emission; negative A is absorption. This deliberately excludes a search over line positions.",
    "For a known line template g, Ahat = sum[g*(flux-1)]/sum(g^2) has standard error s = sigma/sqrt(sum(g^2)). Under H0, A=0 and Ahat ~ Normal(0,s^2). Under H1, the proper prior A ~ Normal(0,w^2) gives Ahat ~ Normal(0,s^2+w^2).",
    "The Bayes factor BF10 is the ratio of those two predictive densities at Ahat. Residual factors orthogonal to the line template cancel. The lower curves are predictive distributions of an estimator, not parameter posteriors.",
    "The reported model probability assumes equal prior model odds. Other prior odds give different posterior model probabilities. The amplitude prior should represent pre-data knowledge, not be tuned after inspecting the line."
  ),
  presets = list("Noise alone" = list(amplitude = 0), "Clear line" = list(amplitude = 2),
                 "Very broad alternative" = list(prior = 10))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
