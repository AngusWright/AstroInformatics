source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "One bad star",
  "Drag a contaminating measurement through an otherwise well-behaved stellar-velocity sample.",
  "A summary's sensitivity to contamination matters, but different summaries can estimate different population properties.",
  controls = list(
    lesson_slider("position", "Contaminant velocity (km/s)", -12, 12, 8, 0.2),
    lesson_slider("contaminants", "Number of replaced measurements", 0, 10, 1, 1),
    lesson_slider("n", "Total sample size", 20, 120, 40, 5)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section1b_html.html",
  section = "1b - Mean, Median, Standard Deviation, and nMAD",
  guide = c("Move one contaminant far from the main population. Which summaries respond?",
            "Increase the number of contaminants while keeping the main population fixed.",
            "Set contamination to zero. Robust does not mean every summary is identical in a finite sample."),
  explanation = c(
    "The uncontaminated population is Normal(0,1) in illustrative velocity units. A chosen number of sample points are replaced by measurements near the contaminant velocity, with SD 0.1. The main Gaussian component is unchanged.",
    "The sample SD uses the n-1 denominator. nMAD is median(abs(x-median(x))) / Phi-inverse(0.75), a Gaussian-consistent robust scale estimate. For a non-Gaussian population, it need not equal its standard deviation.",
    "The horizontal summary bars show mean +/- one sample SD and median +/- one nMAD. They describe location and dispersion, not confidence intervals or uncertainties in the location estimates.",
    "Mean and median answer different questions for an asymmetric population. A contaminant-resistant estimate can target the central population, but ignoring a real secondary population can also discard scientifically important information."
  ),
  presets = list("Clean sample" = list(contaminants = 0),
                 "One extreme star" = list(contaminants = 1, position = 12),
                 "Many contaminants" = list(contaminants = 10, position = 8))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
