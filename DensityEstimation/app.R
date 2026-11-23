source(if (file.exists("teaching.R")) "teaching.R" else "../TeachingTools/teaching.R", local = TRUE)
source("model.R", local = TRUE)

spec <- lesson_spec(
  "Can you trust that bump?",
  "The same stellar-colour sample can look unimodal, bimodal, or full of substructure depending on how it is drawn.",
  "A feature in a density estimate can reflect smoothing or binning rather than a new astrophysical population.",
  controls = list(
    lesson_slider("n", "Number of observed objects", 50, 1000, 200, 50),
    lesson_slider("width", "Histogram bin width", 0.1, 1.5, 0.4, 0.05),
    lesson_slider("origin", "Histogram bin origin (fraction of bin width)", 0, 1, 0.25, 0.05),
    lesson_slider("bandwidth", "Gaussian KDE bandwidth", 0.05, 1, 0.25, 0.05)
  ),
  lecture = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section1b_html.html",
  section = "1b - KDEs and Graphical Distribution Summarisation",
  guide = c("Change bin origin without collecting new data. Do apparent peaks move or disappear?",
            "Compare a very narrow KDE bandwidth with a very broad one.",
            "Increase the sample size while keeping the data-generating population fixed."),
  explanation = c(
    "The known population is a mixture: 60% Normal(-0.8,0.35^2) and 40% Normal(0.9,0.35^2), representing two illustrative stellar-colour populations. The teal curve is this true density, not a fit.",
    "Histogram height is count/(sample size*bin width), so bar area represents probability. Moving the origin changes group membership even though every data value is unchanged.",
    "The KDE places a Gaussian kernel at each observation. Small bandwidths expose sampling fluctuations; large bandwidths can merge real populations. Neither a fixed number of bins nor a single bandwidth is a universal solution.",
    "The displayed KDE peak count ignores peaks below 2% of its maximum height. It is a visual diagnostic, not an estimate with an uncertainty or a formal test for the number of astrophysical components. Kernel tails may extend beyond the finite viewport."
  ),
  presets = list("Undersmoothed" = list(bandwidth = 0.05, width = 0.1),
                 "Oversmoothed" = list(bandwidth = 1, width = 1.5),
                 "A middle ground" = list(bandwidth = 0.25, width = 0.4))
)
parts <- lesson_components(spec, calculate, draw, measure)
ui <- parts$ui
server <- parts$server
shinyApp(ui, server)
