source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Uncertain encounters, many possible galaxies",
  "A statistical companion to the Holmberg simulation: compare a nominal tidal encounter with an ensemble of plausible encounters.",
  "A forward prediction should propagate uncertainty in its inputs, not only show one attractive simulation.",
  controls = list(
    lesson_slider("impact", "Mean impact parameter (model units)", 0.4, 3, 1.5, 0.1),
    lesson_slider("speed", "Mean passage speed (model units)", 0.6, 3, 1.2, 0.1),
    lesson_slider("uncertainty", "Fractional SD of encounter parameters", 0, 0.4, 0.2, 0.02)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2d_html.html",
  section = "2d - Simulation and Monte Carlo Simulations",
  guide = c("Set parameter uncertainty to zero. The ensemble should collapse onto the nominal prediction.",
            "Increase uncertainty and replay the encounter. Which morphological features remain predictable?",
            "Change the mean impact parameter or speed. This changes the physical prediction, not the numerical noise."),
  explanation = c(
    "This is a restricted, collisionless tidal toy model, not a self-consistent N-body galaxy simulation. Forty massless particles orbit a unit-mass softened host potential. A half-mass perturber follows a prescribed straight trajectory. An indirect acceleration keeps the calculation in the host-centred frame.",
    "A kick-drift-kick leapfrog integrates 80 equal time steps of 0.1 in dimensionless units, with softening length 0.3. The perturber is not deflected; self-gravity of the disc, gas, and dynamical friction are absent. All ensemble members are compared at the same physical model time.",
    "Twelve independent encounter parameter pairs are drawn from positive log-normal distributions with the stated arithmetic means and fractional standard deviation. The initial test-particle disc is shared, so the ensemble isolates uncertainty in impact parameter and speed.",
    "The fraction beyond radius 2.5 is a morphological diagnostic, not an escaped or unbound mass fraction. Ensemble spread represents uncertainty under this model and input distribution; it is not automatically a calibrated posterior or a probability for the real Universe.",
    "The original HolmbergSimulation source remains a separate prototype. This page is an astronomical application of the lecture's simulation principles, not an example reproduced from those notes."
  ),
  presets = list("One deterministic encounter" = list(uncertainty = 0),
                 "Uncertain close encounter" = list(impact = 0.7, uncertainty = 0.35)),
  steps = 80
)
parts <- lesson_components(spec, calculate, draw, measure, prepare = prepare)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
