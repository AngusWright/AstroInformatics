source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "From posterior landscape to scientific answer",
  "A narrow, tall peak can contain less probability than a broad, low one.",
  "Marginalisation adds probability across nuisance-parameter values; it does not select the best one.",
  controls = list(
    lesson_slider("weight", "Probability in broad component", 0.05, 0.95, 0.65, 0.05),
    lesson_slider("width", "Broad component standard deviation", 0.5, 2, 1, 0.1),
    lesson_slider("centre", "Centre of highlighted x interval", -2, 3, 1.2, 0.1)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section3b_html.html",
  section = "3b - Multi-parameter Models and Summarising Posterior Information",
  guide = c("Set the broad component to 65% probability. Which peak is taller? Which component contains more probability?",
            "Widen the broad component without changing its probability.",
            "Move the highlighted interval. Its probability is an area under the marginal, not the height of the joint peak."),
  explanation = c(
    "This is a specified, normalised two-component Gaussian mixture used as a toy posterior. The narrow component is centred at (-1,0), with SDs (0.23,0.3); the broad component is centred at (1.3,1), with equal adjustable SDs. Both have correlation 0.45.",
    "The lower curve is the exact marginal p(x) = integral p(x,y) dy over all real y. The highlighted probability integrates that marginal over an x interval of width 0.6. Both calculations are analytic; they do not depend on image resolution.",
    "The upper image has a finite field of view, and its colour uses a square-root stretch to make low-density structure visible. Colour is not a probability per pixel. Contours represent density, not credible-region boundaries.",
    "Profiling maximises over a nuisance parameter; conditioning fixes its value; marginalisation integrates it out. Those are different operations and can produce different answers."
  ),
  presets = list("Tall is not most probable" = list(weight = 0.7, width = 1.5),
                 "Mostly narrow" = list(weight = 0.15, width = 1))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
