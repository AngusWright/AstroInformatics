library(shiny)
library(bslib)

lesson_colours <- c(truth = "#087F8C", data = "#D56534", inference = "#5059A6",
                    muted = "#B7C0CF", text = "#25344B")

lesson_slider <- function(id, label, min, max, value, step = 1) {
  list(id = id, label = label, min = min, max = max, value = value, step = step,
       type = "slider")
}

lesson_choice <- function(id, label, choices, value = unname(choices[1])) {
  list(id = id, label = label, choices = choices, value = value, type = "choice")
}

lesson_spec <- function(title, subtitle, takeaway, controls, lecture, section,
                        guide, explanation, presets = list(), steps = 0) {
  list(title = title, subtitle = subtitle, takeaway = takeaway, controls = controls,
       lecture = lecture, section = section, guide = guide, explanation = explanation,
       presets = presets, steps = steps)
}

lesson_seeded <- function(seed, code) {
  previous <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", .GlobalEnv)
  } else NULL
  on.exit({
    if (is.null(previous)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", previous, .GlobalEnv)
  })
  set.seed(seed)
  force(code)
}

lesson_ribbon <- function(x, lower, upper, colour) {
  polygon(c(x, rev(x)), c(lower, rev(upper)),
          col = adjustcolor(colour, alpha.f = 0.17), border = NA)
}

lesson_histogram <- function(x, breaks) {
  stopifnot(length(breaks) >= 2, all(is.finite(breaks)), all(diff(breaks) > 0),
            all(is.finite(x)))
  index <- findInterval(x, breaks, rightmost.closed = TRUE)
  bins <- length(breaks) - 1
  counts <- tabulate(index[index >= 1 & index <= bins], nbins = bins)
  list(mids = (head(breaks, -1) + tail(breaks, -1)) / 2, counts = counts,
       density = if (length(x)) counts / (length(x) * diff(breaks)) else rep(0, bins),
       excluded = length(x) - sum(counts))
}

lesson_components <- function(spec, calculate, draw, measure, prepare = NULL) {
  controls <- lapply(spec$controls, function(control) {
    if (control$type == "choice") {
      selectInput(control$id, control$label, control$choices, selected = control$value)
    } else {
      sliderInput(control$id, control$label, min = control$min, max = control$max,
                  value = control$value, step = control$step)
    }
  })
  style <- "
    .bslib-page-sidebar {--bslib-sidebar-main-bg: #f5f7fb;}
    .lesson-intro {max-width: 900px; margin: 0 0 1.1rem;}
    .lesson-intro p {font-size: 1.08rem; color: #536176;}
    .lesson-takeaway {background: #e8f3f3; border-left: 4px solid #087f8c;
      border-radius: 0 10px 10px 0; padding: 14px 18px; margin: 1rem 0;}
    .lesson-canvas {background: white; border: 1px solid #e0e5ef;
      border-radius: 14px; padding: 12px; margin: 14px 0;}
    .lesson-metrics {display: flex; flex-wrap: wrap; gap: 12px; margin: 1rem 0;}
    .lesson-metric {flex: 1 1 150px; background: white; border: 1px solid #e0e5ef;
      border-radius: 10px; padding: 12px 16px;}
    .lesson-metric strong {display: block; font-size: 1.3rem; color: #25344b;}
    .lesson-metric span {font-size: .86rem; color: #536176;}
    .lesson-source {border-top: 1px solid #dfe5ef; padding: 16px 0; margin-top: 20px;}
    .lesson-controls {display: flex; gap: 6px; flex-wrap: wrap; margin: 12px 0;}
    .lesson-preset {margin: 3px 4px 3px 0;}
    .lesson-explain {max-width: 850px; line-height: 1.65; padding-top: 1rem;}
    .nav-tabs {margin-top: 10px;}
    @media (max-width: 700px) {.lesson-canvas {padding: 2px;}}
  "
  ui <- page_sidebar(
    title = spec$title, fillable = FALSE,
    theme = bs_theme(version = 5, bg = "#F5F7FB", fg = "#25344B", primary = "#405B99"),
    sidebar = sidebar(
      open = "open", width = 300,
      h5("Explore"), controls,
      if (length(spec$presets)) tagList(
        h6("Guided examples"),
        lapply(seq_along(spec$presets), function(i) {
          actionButton(paste0("preset_", i), names(spec$presets)[i], class = "lesson-preset")
        })
      ),
      hr(), actionButton("redraw", "New realisation"),
      tags$details(
        tags$summary("Simulation settings"),
        numericInput("seed", "Random seed", 42, min = 0, max = 10000000, step = 1),
        if (spec$steps > 0) selectInput("delay", "Animation speed",
                                       c("Fast" = 100, "Normal" = 200, "Slow" = 500),
                                       selected = 200),
        helpText("Controls reuse the same underlying random draws. New realisation changes the draws; animation only reveals them.")
      ),
      hr(), a("All demonstrations", href = "../"),
      p(a("Read the lecture", href = spec$lecture))
    ),
    tags$head(tags$style(HTML(style))),
    div(class = "lesson-intro", p(spec$subtitle)),
    div(class = "lesson-takeaway", strong("The principle: "), spec$takeaway),
    tabsetPanel(
      tabPanel(
        "Explore",
        if (spec$steps > 0) tagList(
          div(class = "lesson-controls",
              actionButton("play", "Play animation", class = "btn-primary"),
              actionButton("step", "Step"), actionButton("rewind", "Rewind"),
              actionButton("finish", "Show result")),
          uiOutput("progress")
        ),
        div(class = "lesson-canvas", plotOutput("plot", height = "500px")),
        uiOutput("metrics"),
        p(class = "text-muted", "A controlled, simulated teaching example; not a fit to a real survey.")
      ),
      tabPanel(
        "Explain",
        div(class = "lesson-explain",
            h4("Try this"),
            tags$ol(lapply(spec$guide, tags$li)),
            h4("What the demonstration assumes"),
            lapply(spec$explanation, p))
      )
    ),
    div(class = "lesson-source",
        strong("Read the lecture: "),
        a(spec$section, href = spec$lecture),
        span("  |  "), a("All demonstrations", href = "../"))
  )

  server <- function(input, output, session) {
    frame <- reactiveVal(spec$steps)
    playing <- reactiveVal(FALSE)
    generation <- reactiveVal(0L)
    controls_changed <- reactive({
      c(lapply(spec$controls, function(control) input[[control$id]]),
        list(input$seed, generation()))
    })
    observeEvent(controls_changed(), {
      playing(FALSE)
      frame(spec$steps)
    }, priority = 100)
    observeEvent(input$redraw, generation(generation() + 1L))
    for (i in seq_along(spec$presets)) local({
      index <- i
      observeEvent(input[[paste0("preset_", index)]], {
        values <- spec$presets[[index]]
        for (control in spec$controls) if (control$id %in% names(values)) {
          if (control$type == "choice") {
            updateSelectInput(session, control$id, selected = values[[control$id]])
          } else {
            updateSliderInput(session, control$id, value = values[[control$id]])
          }
        }
      })
    })
    observeEvent(input$play, {
      if (playing()) {
        playing(FALSE)
      } else {
        if (frame() >= spec$steps) frame(0L)
        playing(TRUE)
      }
    })
    observeEvent(input$step, {
      playing(FALSE)
      frame(min(frame() + 1L, spec$steps))
    })
    observeEvent(input$rewind, {
      playing(FALSE)
      frame(0L)
    })
    observeEvent(input$finish, {
      playing(FALSE)
      frame(spec$steps)
    })
    observe({
      updateActionButton(session, "play", label = if (playing()) "Pause" else "Play animation")
    })
    observe({
      if (!playing()) return()
      invalidateLater(as.numeric(input$delay), session)
      isolate({
        frame(min(frame() + 1L, spec$steps))
        if (frame() == spec$steps) playing(FALSE)
      })
    })
    check_seed <- function() {
      req(all(vapply(spec$controls, function(control) !is.null(input[[control$id]]), logical(1))))
      validate(need(length(input$seed) == 1 && is.finite(input$seed) &&
                      input$seed >= 0 && input$seed <= 10000000 &&
                      input$seed == round(input$seed), "Enter an integer seed from 0 to 10000000."))
    }
    prepared <- reactive({
      controls_changed()
      check_seed()
      lesson_seeded(input$seed + generation(), prepare(input))
    })
    result <- reactive({
      check_seed()
      phase <- if (spec$steps > 0) frame() / spec$steps else 1
      if (is.null(prepare)) {
        lesson_seeded(input$seed + generation(), calculate(input, phase))
      } else {
        calculate(input, phase, prepared())
      }
    })
    output$plot <- renderPlot({
      par(mar = c(4.5, 4.5, 2.5, 1.2), bg = "white", fg = lesson_colours["text"],
          col.axis = "#536176", col.lab = "#25344B", family = "sans",
          las = 1, bty = "l", cex.main = 1.1)
      draw(result(), input)
    }, res = 110, alt = paste(spec$title, "-", spec$takeaway))
    output$metrics <- renderUI({
      values <- measure(result(), input)
      div(class = "lesson-metrics", lapply(seq_along(values), function(i) {
        div(class = "lesson-metric", strong(unname(values[i])), span(names(values)[i]))
      }))
    })
    output$progress <- renderUI({
      tagList(
        tags$progress(value = frame(), max = spec$steps, style = "width:100%;",
                      `aria-label` = "Illustration progress"),
        span(class = "text-muted", sprintf(" Illustrated step %d of %d. %s",
                                           frame(), spec$steps,
                                           if (playing()) "Playing." else "Paused; Play reveals the same realisation in order."))
      )
    })
  }
  list(ui = ui, server = server)
}
