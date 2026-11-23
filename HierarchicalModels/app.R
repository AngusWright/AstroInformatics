source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Many noisy stars, one population",
  "Infer individual stellar velocities and a shared population mean without pretending that every star is identical.",
  "Partial pooling borrows population information most strongly for the least precisely measured objects.",
  controls = list(
    lesson_choice("pooling", "Inference model", c("Partial pooling" = "partial", "No pooling" = "none", "Complete pooling" = "complete")),
    lesson_slider("scatter", "Intrinsic velocity dispersion (km/s)", 0, 3, 1, 0.1),
    lesson_slider("noise", "Measurement uncertainty multiplier", 0.3, 3, 1, 0.1)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section4c_html.html",
  section = "4c - Hierarchical Models and Building a Model",
  guide = c("Compare no pooling and partial pooling. Which measurements move most?",
            "Increase observational noise while holding the population dispersion fixed.",
            "Compare complete pooling with a population that really has intrinsic scatter. A simpler model can be too restrictive."),
  explanation = c(
    "Fourteen member-star velocities are drawn from a Normal population with true mean 1.5 km/s. Each star has its own known measurement standard deviation. Membership is known and selection is ignored.",
    "The partial-pooling model is theta_i | mu ~ Normal(mu, tau^2), x_i | theta_i ~ Normal(theta_i, s_i^2), and mu ~ Normal(0, 5^2). The dispersion tau is known in this demonstration. The population mean mu is inferred, and its uncertainty is propagated into every individual interval.",
    "For V = 1/(1/25 + sum[1/(tau^2+s_i^2)]), m = V sum[x_i/(tau^2+s_i^2)], and w_i = tau^2/(tau^2+s_i^2), the individual posterior mean is w_i*x_i + (1-w_i)*m. Its variance is tau^2*s_i^2/(tau^2+s_i^2) + (1-w_i)^2*V.",
    "Complete pooling fits a separate tau=0 model. No pooling uses separate flat individual priors, giving Normal(x_i,s_i^2) posteriors. During animation, positions and interval widths interpolate for illustration; only the final intervals are the stated posterior intervals."
  ),
  presets = list("Precise stars" = list(noise = 0.3, pooling = "partial"),
                 "Uncertain stars" = list(noise = 2.5, pooling = "partial")),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
