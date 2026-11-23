source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "The universe you don't detect",
  "A telescope's depth can create a luminosity-distance relation in a population with no intrinsic evolution.",
  "An observed catalogue is not the parent population. Selection can imitate a physical trend.",
  controls = list(
    lesson_slider("limit", "Survey limiting apparent magnitude", 17, 24, 20, 0.1),
    lesson_choice("overlay", "Reveal the missing population", c("Show" = "show", "Hide" = "hide"))
  ),
  lecture = "https://anguswright.github.io/LectureNotes/AstroStatistics_1.html",
  section = "Astrostatistics 1 - Detection Effects, Malmquist Bias, and Survivor Bias",
  guide = c("Use the shallow-survey preset. At large distances, which galaxies survive?",
            "Reveal the missing population: was the underlying relation ever present?",
            "Replay the animation to watch the catalogue change as the flux limit becomes more restrictive."),
  explanation = c(
    "This Euclidean toy universe has a distance-independent log-normal luminosity population and uniform spatial density. Distance is in Mpc; luminosity is in arbitrary reference units. The magnitude relation is m = M + 5 log10(distance/Mpc) + 25, with M = -18 - 2.5 log10(L/L0). Cosmological evolution and K-corrections are deliberately excluded.",
    "Objects are detected if m is below the survey limit. The orange curve is the sample mean log-luminosity in distance bins; the dashed teal line is the true population expectation, zero. Empty bins have no estimated mean.",
    "The appropriate selected-data model is p(x | detected) = S(x) p(x) / integral[S(u) p(u) du]. The denominator depends on the population model; simply fitting the detected points can therefore be misleading."
  ),
  presets = list("Shallow survey" = list(limit = 18.5), "Deep survey" = list(limit = 23)),
  steps = 60
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
