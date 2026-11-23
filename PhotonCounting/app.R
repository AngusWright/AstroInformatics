source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Counting photons",
  "Individual photons arrive unpredictably, but their count distribution follows a simple law.",
  "Poisson fluctuations are discrete and asymmetric at low counts; a Gaussian is only an approximation.",
  controls = list(
    lesson_slider("rate", "Photon arrival rate (per second)", 0.2, 8, 2, 0.2),
    lesson_slider("exposure", "Exposure time (seconds)", 0.5, 20, 5, 0.5),
    lesson_choice("gaussian", "Gaussian approximation", c("Show" = "show", "Hide" = "hide"))
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2c_html.html",
  section = "2c - The Poisson Distribution and Gaussian Distribution",
  guide = c("Choose a faint source and a short exposure. Can a Gaussian assign probability to impossible counts?",
            "Replay to accumulate photons over time.",
            "Increase exposure and compare the count SD with the expected count and relative precision."),
  explanation = c(
    "The top panel is one homogeneous Poisson arrival process with independent exponential waiting times. There is no background, detector dead time, or varying source rate.",
    "At exposure t and known rate r, N ~ Poisson(lambda=r*t), with expectation and variance both lambda. The relative count SD is 1/sqrt(lambda) for positive lambda. This is a sampling model, not an uncertainty interval inferred from one realised count.",
    "The lower histogram represents 3,000 independent Poisson counts at the current exposure. Fixed uniform quantiles couple successive histogram frames monotonically; these histogram quantiles are not 3,000 simulated event timelines.",
    "The Gaussian approximation has mean lambda and SD sqrt(lambda). Its probabilities are integrated over each unit-width integer bin (continuity correction), so both curves show comparable probability masses. Its area below -0.5 corresponds to impossible negative integer counts after rounding. At zero exposure the distribution is a point mass at zero, and no Gaussian approximation or relative precision is defined."
  ),
  presets = list("Few photons" = list(rate = 0.4, exposure = 1),
                 "Many photons" = list(rate = 8, exposure = 15)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
