source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Why stacking works - and when it doesn't",
  "Build distributions of averages from skewed, mixed, or extremely heavy-tailed observations.",
  "The central limit theorem concerns the sampling distribution of an average, not the shape of the original data.",
  controls = list(
    lesson_choice("population", "Observation/error population",
                  c("Skewed exponential" = "exponential", "Rare bright outliers" = "mixture", "Cauchy: no finite variance" = "cauchy")),
    lesson_slider("n", "Measurements in each average", 1, 100, 30, 1),
    lesson_slider("repetitions", "Independent averaging experiments", 200, 3000, 1000, 100)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2c_html.html",
  section = "2c - Central Limit Theorem and the Student t Distribution",
  guide = c("Start with one exponential observation per average, then replay to grow each stack.",
            "Try a population with occasional bright outliers. More repetitions improve the displayed histogram; more measurements change the sampling distribution.",
            "Select Cauchy errors. A bigger stack does not create a finite mean or variance."),
  explanation = c(
    "Each experiment averages n independent identically distributed measurements. Exponential observations have mean 1 and variance 1. The mixture is 95% Normal(0,1) plus 5% Normal(5,1), with mean 0.25 and variance 2.1875.",
    "For finite variance, the lower panel plots sqrt(n)*(sample mean - population mean)/population SD. The ordinary iid CLT predicts convergence to Normal(0,1). This standardisation is essential: the raw mean distribution gets narrower rather than retaining unit variance.",
    "Cauchy observations are Student-t with one degree of freedom. Neither their mean nor variance exists; their sample mean has the same Cauchy distribution for every n. The lower panel therefore shows unstandardised Cauchy means and an exact Cauchy reference, not a Gaussian.",
    "The finite plotting windows omit extreme tails without renormalising the histogram. Stacking does not remove shared calibration errors, correlations, selection bias, or unmodelled variability. The observations here deliberately exclude those effects."
  ),
  presets = list("One measurement" = list(n = 1), "Large stack" = list(n = 100),
                 "A counterexample" = list(population = "cauchy", n = 100)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure, prepare = prepare)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
