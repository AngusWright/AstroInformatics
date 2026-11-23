source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "A sky full of false discoveries",
  "Every resolution element here contains only noise. Some will still look like sources.",
  "A rare result at one chosen location can be common somewhere in a large search.",
  controls = list(
    lesson_slider("side", "Independent resolution elements per side", 5, 60, 40, 1),
    lesson_slider("threshold", "Detection threshold (Gaussian sigma)", 1.5, 5, 3, 0.1),
    lesson_choice("tails", "What counts as a detection?", c("Positive peaks" = "one", "Either sign" = "two"))
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section4a_html.html",
  section = "4a - Significance, Variable Selection, and the Problem of Choice",
  guide = c("Keep the threshold fixed and increase the area searched.",
            "Replay the search: no real source is ever injected.",
            "Compare the single-location false-positive probability with the probability of at least one false positive anywhere."),
  explanation = c(
    "Each displayed cell represents one independent standard-Normal resolution element, not necessarily one camera pixel. Real PSFs and resampling correlate pixels, so real trial counts require more care than counting pixels.",
    "For a positive threshold z, the single-location probability is p = 1-Phi(z); a two-sided search doubles it. For N independent locations, P(at least one false positive) = 1-(1-p)^N. The expected number of false positives is N*p, which is not that probability.",
    "A p-value is not the probability that a source is false. This page calibrates a specified noise-only search procedure. Choosing thresholds, masks, or additional searches after inspecting the image changes the procedure."
  ),
  presets = list("Small search" = list(side = 5, threshold = 3),
                 "Large search" = list(side = 60, threshold = 3)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
