source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Correlations that disappear",
  "Two independent intrinsic galaxy properties can look related when a shared distance effect is left in the data.",
  "An observed association does not establish a direct causal relationship.",
  controls = list(
    lesson_slider("scatter", "Independent intrinsic scatter in each band", 0.1, 2, 0.4, 0.1),
    lesson_slider("centre", "Centre of normalised log-distance slice", 0, 1, 0.5, 0.05),
    lesson_slider("width", "Width of distance slice", 0.05, 1, 0.2, 0.05)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/AstroStatistics_1.html",
  section = "Astrostatistics 1 - Interpreting Correlation and Confounding Variables",
  guide = c("Start with the whole distance range. Both bands share a distance dependence.",
            "Narrow the distance slice, highlighting fewer objects. Compare the full-sample and slice correlations.",
            "Increase intrinsic scatter. Does a common cause always produce an obvious correlation?"),
  explanation = c(
    "The toy log-fluxes are X=-2Z+epsilon_x and Y=-2Z+epsilon_y, with Z uniform on [0,1] and independent Normal intrinsic terms of known SD. Z represents a rescaled log-distance; the flux zero points are arbitrary. There is no direct X-to-Y link in this specified model.",
    "Conditioning on an exact distance makes X and Y independent here. A finite-width distance slice still contains distance variation, so its sample correlation need not be zero. Sampling fluctuations also remain.",
    "This demonstration establishes a common cause because its generating model is known. In real observations, conditioning on a convenient variable does not by itself establish causality and can introduce collider bias. Flux limits and measurement errors are intentionally absent."
  ),
  presets = list("Whole population" = list(width = 1, centre = 0.5),
                 "Narrow distance slice" = list(width = 0.05, centre = 0.5))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
