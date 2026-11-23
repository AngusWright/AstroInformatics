source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Turning uniform randomness into a universe",
  "Map equally likely uniform draws through a quantile function to assemble a non-uniform population.",
  "If U is Uniform(0,1), then F-inverse(U) follows the target distribution.",
  controls = list(
    lesson_choice("family", "Positive target distribution",
                  c("Truncated stellar-mass power law" = "power", "Exponential" = "exponential", "Log-normal" = "lognormal")),
    lesson_slider("slope", "Power-law exponent alpha (power-law family only)", 0.5, 3.5, 2.35, 0.05),
    lesson_slider("n", "Number of random draws", 50, 1500, 500, 50)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2d_html.html",
  section = "2d - Inversion Sampling",
  guide = c("Replay and follow the horizontal uniform-probability guides to the CDF, then down to the sampled value.",
            "Change the power-law exponent. Uniform draws remain uniform; transformed draws do not.",
            "Increase the number of draws. The histogram becomes smoother without changing the target density."),
  explanation = c(
    "The power-law density is proportional to x^(-alpha) on [0.1,20], in illustrative solar-mass units. Alpha=1 uses a logarithmic CDF and inverse. The exponential has rate 1; the log-normal has log-mean 0 and log-SD 0.7.",
    "All distributions are positive and displayed in log10(x). The CDF remains F(10^z), but the lower density must include the Jacobian: p_z(z) = ln(10)*10^z*p_x(10^z). A density on one scale is not the same curve on a transformed scale.",
    "Thin guides show up to twelve of the most recently revealed draws. Histogram probabilities are normalised by all revealed draws. With zero draws, only the analytic target is shown.",
    "The reported discrepancy is the empirical CDF's maximum absolute difference from the target (the one-sample KS statistic), not a significance or model-selection probability."
  ),
  presets = list("Log-uniform power law" = list(family = "power", slope = 1),
                 "Steep stellar population" = list(family = "power", slope = 3)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
