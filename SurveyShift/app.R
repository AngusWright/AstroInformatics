source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "A classifier under alien skies",
  "A classifier that works on one survey can fail when another survey changes the measured colours.",
  "Validation is only informative for deployment conditions that it actually represents.",
  controls = list(
    lesson_slider("offset", "New-survey colour calibration offset", -2, 2, 0.8, 0.1),
    lesson_slider("noise", "Additional new-survey colour noise", 0, 1.2, 0.3, 0.05),
    lesson_choice("view", "Display validation sample", c("New survey" = "target", "Original survey" = "original"))
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section4c_html.html",
  section = "4c - Machine Learning, Extrapolation, and Covariate Shift",
  guide = c("Set offset and additional noise to zero: the two surveys follow the same model.",
            "Change the new survey's calibration without retraining the decision rule.",
            "Switch the displayed sample. Compare independent original-survey validation with new-survey performance."),
  explanation = c(
    "The toy training sample has 200 stars and 200 quasars in a two-colour diagram. A linear discriminant is fitted using sample class means and a pooled within-class covariance matrix. No external machine-learning package is used.",
    "Both validation samples are independent of training and contain equally many objects from each class. Target observations add a calibration offset (delta,0.4*delta) and independent Gaussian colour noise. The astrophysical class labels do not change.",
    "This is general survey/domain shift. Changing the measurement process can change P(class | measured colour), so this is not restricted to the narrower mathematical definition of pure covariate shift.",
    "The background shows the frozen classifier's prediction, not a probability of correctness. Orange rings reveal mistakes using simulated truth. Neither balanced class frequencies nor these simple colour populations are intended as a realistic quasar survey."
  ),
  presets = list("Matched surveys" = list(offset = 0, noise = 0),
                 "Calibration mismatch" = list(offset = 1.5, noise = 0.3))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
