source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Why reliable detectors still find impostors",
  "A rare-event survey can produce more false candidates than genuine discoveries, even with a good detector.",
  "P(real event | flagged) is not P(flagged | real event). The population base rate matters.",
  controls = list(
    lesson_slider("log_prevalence", "log10 fraction of real events", -4, -1, -2, 0.25),
    lesson_slider("sensitivity", "Detection efficiency for real events", 0.5, 1, 0.95, 0.01),
    lesson_slider("false_positive", "False-positive rate for non-events", 0.001, 0.1, 0.01, 0.001),
    lesson_choice("reveal", "Reveal simulated candidate identities", c("Reveal truth" = "show", "Hide truth" = "hide"))
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2b_html.html",
  section = "2b - Conditional Probability, The Prosecutor's Fallacy, and Anomaly Detection",
  guide = c("Keep the detector fixed and make genuine events rarer.",
            "Compare sensitivity with the fraction of flagged candidates that are genuine.",
            "Hide the simulated truth: a real detector only supplies candidates, not their true identities."),
  explanation = c(
    "The simulated survey has 10,000 independent objects. A real event has probability q; its flag probability is sensitivity s. A non-event has false-positive probability f. These are known model inputs, not estimated performance measurements.",
    "Bayes' rule gives P(real | flag) = q*s / [q*s + (1-q)*f]. Sensitivity is P(flag | real), a different conditional probability. The false-positive rate is not the fraction of all flags that are false.",
    "Candidate identities are visible only because this is a simulation. The displayed true fraction fluctuates around its theoretical value, especially for rare events and small candidate samples. No flags is a legitimate outcome, not a zero posterior probability."
  ),
  presets = list("Rare events" = list(log_prevalence = -4),
                 "Commoner events" = list(log_prevalence = -1)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
