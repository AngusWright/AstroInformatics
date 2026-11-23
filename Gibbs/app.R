library(shiny)
library(bslib)
source("gibbs.R", local = TRUE)

ui <- page_sidebar(
  title = "Gibbs sampling: one conditional at a time",
  fillable = FALSE,
  sidebar = sidebar(
    open = "open", width = 330,
    p(a("Read the lecture: The Gibbs Sampler",
        href = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section3c_html.html#building-an-mcmc-sampler-the-gibbs-sampler")),
    p(a("All demonstrations", href = "../")),
    accordion(
      accordion_panel(
        "Likelihood and target",
        selectInput("family", "Likelihood family",
                    c("Multivariate Gaussian (2D)" = "gaussian",
                      "Multivariate Student-t (2D)" = "student",
                      "Poisson with Gamma hierarchy" = "poisson")),
        conditionalPanel(
          "input.family != 'poisson'",
          sliderInput("location1", "Observation y1 (posterior location 1)",
                      -3, 3, 0, step = 0.1),
          sliderInput("location2", "Observation y2 (posterior location 2)",
                      -3, 3, 0, step = 0.1),
          sliderInput("scale1", "Scale 1 (Gaussian SD; Student-t scale)",
                      0.25, 3, 1, step = 0.05),
          sliderInput("scale2", "Scale 2 (Gaussian SD; Student-t scale)",
                      0.25, 3, 1, step = 0.05),
          sliderInput("rho", "Correlation rho", -0.98, 0.98, 0.8, step = 0.01)
        ),
        conditionalPanel(
          "input.family == 'student'",
          sliderInput("df", "Student-t degrees of freedom (nu)",
                      3, 30, 5, step = 1)
        ),
        conditionalPanel(
          "input.family == 'poisson'",
          sliderInput("n", "Number of unit-exposure observations", 1, 100, 10, step = 1),
          sliderInput("total", "Observed total count", 0, 200, 30, step = 1),
          sliderInput("a", "Shape a of lambda | beta", 1, 10, 2, step = 0.5),
          sliderInput("hyper_shape", "Hyperprior Gamma shape c", 1, 10, 2, step = 0.5),
          sliderInput("hyper_rate", "Hyperprior Gamma rate d", 0.1, 5, 1, step = 0.1)
        )
      ),
      accordion_panel(
        "Chain settings",
        sliderInput("sweeps", "Number of complete sweeps", 10, 2000, 100, step = 10),
        numericInput("seed", "Random seed (change for a new trajectory)", 666,
                     min = 0, max = .Machine$integer.max, step = 1),
        conditionalPanel(
          "input.family != 'poisson'",
          numericInput("start1", "Starting theta1", -2.5, step = 0.5),
          numericInput("start2", "Starting theta2", 2.5, step = 0.5)
        ),
        conditionalPanel(
          "input.family == 'poisson'",
          numericInput("start_lambda", "Starting lambda (> 0)", 1, min = 0.01, step = 0.1),
          numericInput("start_beta", "Starting beta (> 0)", 1, min = 0.01, step = 0.1)
        ),
        helpText("Changing the target, seed, starting point, or sweep count generates a new chain and rewinds playback.")
      ),
      accordion_panel(
        "Playback and display",
        selectInput("delay", "Time per coordinate update",
                    c("Fast (0.1 s)" = 100, "0.25 s" = 250,
                      "0.5 s" = 500, "Slow (1 s)" = 1000), selected = 500),
        sliderInput("trail", "Recent coordinate moves to show", 10, 200, 40, step = 10),
        sliderInput("burn", "Discard initial complete sweeps (burn-in)",
                    0, 99, 10, step = 1),
        checkboxInput("show_samples", "Show completed-sweep samples", TRUE),
        helpText("Playback, trail length, and burn-in never change the generated trajectory.")
      ),
      open = c("Likelihood and target", "Chain settings")
    )
  ),
  p("Gibbs sampling builds a joint sample by drawing one parameter from its full conditional, ",
    "holding the others fixed. Here two coordinate updates make one complete sweep."),
  tabsetPanel(
    tabPanel(
      "Animated updates",
      div(
        actionButton("play", "Play", class = "btn-primary"),
        actionButton("step", "Next coordinate"),
        actionButton("sweep", "Finish / next sweep"),
        actionButton("rewind", "Rewind"),
        actionButton("finish", "Show all iterations")
      ),
      fluidRow(
        column(4, numericInput("jump_to", "Jump to coordinate update", 0, min = 0, max = 200, step = 1)),
        column(3, br(), actionButton("jump", "Jump"))
      ),
      uiOutput("progress"),
      uiOutput("step_description"),
      fluidRow(
        column(8, plotOutput("joint_plot", height = "470px")),
        column(4, plotOutput("conditional_plot", height = "350px"),
               uiOutput("conditional_description"))
      ),
      p("Grey contours show relative joint posterior density, not credible regions. ",
        "Orange horizontal moves update parameter 1; blue vertical moves update parameter 2. ",
        "Open circles mark the previous state; the filled circle is the current state. ",
        "Small green dots are retained complete sweeps; grey dots are burn-in."),
      tableOutput("step_table")
    ),
    tabPanel(
      "Trace and marginals",
      p("Only fully completed sweeps up to the displayed frame enter these summaries. ",
        "Intermediate coordinate states and the initial state are not counted as separate retained sweeps."),
      uiOutput("sample_count"),
      plotOutput("diagnostic_plot", height = "650px"),
      tableOutput("sample_summary"),
      p("Solid marginal curves are the exact Gaussian or Student-t posterior marginals. ",
        "For the Poisson hierarchy, histograms are shown without an analytic marginal overlay. ",
        "Finite, correlated samples need not match the target closely."),
      plotOutput("acf_plot", height = "310px"),
      p("High autocorrelation means many sweeps may carry little new information. ",
        "Burn-in discards early iterations; it does not remove dependence or certify convergence.")
    ),
    tabPanel(
      "How Gibbs works",
      h4("A systematic-scan Gibbs sweep"),
      tags$ol(
        tags$li("Choose an initial state (theta1, theta2) within the target support."),
        tags$li("Draw theta1 from p(theta1 | current theta2, data). Hold theta2 fixed: a horizontal move."),
        tags$li("Draw theta2 from p(theta2 | NEW theta1, data). Hold the newly sampled theta1 fixed: a vertical move."),
        tags$li("Record the completed pair, then repeat. Each sweep uses the most recently available coordinate values.")
      ),
      tags$pre("theta1[t+1] ~ p(theta1 | theta2[t], data)\ntheta2[t+1] ~ p(theta2 | theta1[t+1], data)"),
      p("In d dimensions, update each coordinate from its full conditional given all the others, ",
        "using the newly updated coordinates immediately. This page visualises the two-dimensional case."),
      h4("Why it works"),
      p("Each exact full-conditional update preserves the joint target distribution. ",
        "Under suitable irreducibility and convergence conditions, repeated sweeps approach that target. ",
        "Gibbs draws are not optimisation steps: moves to lower joint density are allowed. ",
        "There is no Metropolis accept/reject step when sampling these full conditionals exactly."),
      h4("Models and full conditionals"),
      uiOutput("model_equations"),
      h4("Experiments"),
      tags$ol(
        tags$li("Start with rho = 0.8, y = (0, 0), scales = 1, and a dispersed starting point, as in the lecture. ",
                "Step through the first sweep and watch the conditioning value change."),
        tags$li("Compare rho = 0 with rho = 0.98. Strong correlation makes the conditional narrow, ",
                "so the chain moves slowly along the long axis of the posterior. Try negative correlation too."),
        tags$li("Select Student-t. Moving far from the centre widens its full conditional, unlike the Gaussian. ",
                "The conditional has nu + 1 degrees of freedom, not nu."),
        tags$li("Select the Poisson hierarchy. Both coordinates must remain positive, and each draw is Gamma. ",
                "Changing the count data changes the posterior, not just its visual appearance."),
        tags$li("Use 'Show all iterations' to examine traces and marginals. Change the starting point and seed ",
                "to compare independent runs: one attractive trace or histogram is not proof of convergence.")
      ),
      p("The seed fixes the complete trajectory. Animation reveals those draws in order; ",
        "pausing, replaying, or changing the display does not draw a new sample. Playback starts paused."),
      a("Lecture notes: Building an MCMC Sampler: The Gibbs Sampler",
        href = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section3c_html.html#building-an-mcmc-sampler-the-gibbs-sampler")
    )
  )
)

server <- function(input, output, session) {
  target <- reactive({
    req(input$family)
    gibbs_target(input$family, location = c(input$location1, input$location2),
                 scale = c(input$scale1, input$scale2), rho = input$rho, df = input$df,
                 n = input$n, total = input$total, a = input$a,
                 c = input$hyper_shape, d = input$hyper_rate)
  })
  chain <- reactive({
    model <- target()
    start <- if (model$positive) c(input$start_lambda, input$start_beta) else c(input$start1, input$start2)
    validate(
      need(length(start) == 2 && all(is.finite(start)), "Enter two finite starting coordinates."),
      need(!model$positive || all(start > 0), "Poisson starting coordinates must both be positive."),
      need(length(input$seed) == 1 && is.finite(input$seed) &&
             input$seed >= 0 && input$seed <= .Machine$integer.max &&
             input$seed == round(input$seed), "Enter an integer seed between 0 and 2147483647.")
    )
    gibbs_chain(model, start, input$sweeps, input$seed)
  })
  frame <- reactiveVal(0L)
  playing <- reactiveVal(FALSE)
  observeEvent(chain(), {
    playing(FALSE)
    frame(0L)
    updateNumericInput(session, "jump_to", value = 0, max = 2 * chain()$sweeps)
    updateSliderInput(session, "burn", max = chain()$sweeps - 1,
                      value = min(isolate(input$burn), chain()$sweeps - 1))
  }, priority = 100)

  seek <- function(update) {
    playing(FALSE)
    frame(as.integer(update))
  }
  observeEvent(input$play, {
    if (playing()) {
      playing(FALSE)
    } else {
      if (frame() == 2 * chain()$sweeps) frame(0L)
      playing(TRUE)
    }
  })
  observeEvent(input$step, seek(min(frame() + 1L, 2 * chain()$sweeps)))
  observeEvent(input$sweep, seek(min(2L * (floor(frame() / 2) + 1L), 2 * chain()$sweeps)))
  observeEvent(input$rewind, seek(0L))
  observeEvent(input$finish, seek(2 * chain()$sweeps))
  observeEvent(input$jump, {
    value <- input$jump_to
    if (length(value) != 1 || !is.finite(value) || value != round(value) ||
        value < 0 || value > 2 * chain()$sweeps) {
      showNotification(sprintf("Choose an integer update from 0 to %d.", 2 * chain()$sweeps), type = "error")
    } else {
      seek(value)
    }
  })
  observe({
    updateActionButton(session, "play", label = if (playing()) "Pause" else "Play")
  })
  observe({
    if (!playing()) return()
    invalidateLater(as.numeric(input$delay), session)
    isolate({
      next_frame <- min(frame() + 1L, 2 * chain()$sweeps)
      frame(next_frame)
      if (next_frame == 2 * chain()$sweeps) playing(FALSE)
    })
  })

  current <- reactive(gibbs_frame(chain(), frame()))
  completed <- reactive(gibbs_samples(chain(), frame()))
  retained <- reactive(gibbs_samples(chain(), frame(), input$burn))
  colours <- c("#D55E00", "#0072B2")

  grid <- reactive({
    model <- chain()$target
    limits <- lapply(1:2, function(j) {
      if (model$positive) c(0, 1.05 * max(chain()$states[, j])) else range(model$limits[[j]], chain()$states[, j])
    })
    x <- seq(limits[[1]][1], limits[[1]][2], length.out = 180)
    y <- seq(limits[[2]][1], limits[[2]][2], length.out = 180)
    if (model$positive) {
      # Resolve a narrow posterior even when the starting point is far away.
      probabilities <- seq(0.001, 0.999, length.out = 120)
      x <- sort(unique(c(x[x > 0], quantile(chain()$states[, 1], probabilities))))
      y <- sort(unique(c(y[y > 0], quantile(chain()$states[, 2], probabilities))))
    }
    z <- outer(x, y, model$log_density)
    list(x = x, y = y, z = z, limits = limits,
         levels = max(z) + log(c(0.01, 0.05, 0.2, 0.5, 0.8)))
  })

  output$progress <- renderUI({
    total <- 2 * chain()$sweeps
    tagList(
      p(strong(sprintf("Update %d / %d | Completed sweeps: %d / %d",
                       frame(), total, floor(frame() / 2), chain()$sweeps))),
      tags$progress(value = frame(), max = total, style = "width: 100%;",
                    `aria-label` = "Coordinate updates completed")
    )
  })
  output$step_description <- renderUI({
    view <- current()
    other <- 3L - view$coordinate
    p(strong(if (view$has_draw) sprintf("Last draw: parameter %d.", view$coordinate) else "Ready: draw parameter 1."),
      sprintf(" Hold %s fixed at %.4g. %s",
              chain()$target$axes[other], view$before[other],
              if (view$has_draw) sprintf("New value: %.4g.", view$after[view$coordinate]) else ""),
      if (frame() < 2 * chain()$sweeps) {
        sprintf(" Next: parameter %d, conditioned on the current other coordinate.", frame() %% 2L + 1L)
      } else " All sweeps complete; Play replays this same chain.")
  })
  output$conditional_description <- renderUI({
    tagList(
      h5(if (current()$has_draw) "Full conditional used for the last draw" else "Full conditional for the first draw"),
      p(current()$distribution$label),
      p("This is a distribution for one coordinate with the other fixed, not a marginal distribution.")
    )
  })
  output$joint_plot <- renderPlot({
    model <- chain()$target
    background <- grid()
    view <- current()
    par(mar = c(4, 4, 2, 1))
    plot(NA, xlim = background$limits[[1]], ylim = background$limits[[2]],
         xlab = model$axes[1], ylab = model$axes[2], main = model$label,
         asp = if (model$positive) NA_real_ else 1)
    contour(background$x, background$y, background$z, levels = background$levels,
            add = TRUE, drawlabels = FALSE, col = "#BBBBBB")
    if (isTRUE(input$show_samples) && nrow(completed()) > 0) {
      points(completed(), pch = 16, cex = 0.6,
             col = ifelse(seq_len(nrow(completed())) <= input$burn, "#99999988", "#009E7366"))
    }
    start <- max(0, frame() - input$trail)
    path <- chain()$states[seq.int(start + 1L, frame() + 1L), , drop = FALSE]
    lines(path, col = "#55555599", lwd = 1.3)
    points(chain()$states[1, 1], chain()$states[1, 2], pch = 4, cex = 1.5, lwd = 2)
    colour <- colours[view$coordinate]
    if (view$coordinate == 1L) abline(h = view$before[2], col = colour, lty = 3) else abline(v = view$before[1], col = colour, lty = 3)
    if (view$has_draw) {
      arrows(view$before[1], view$before[2], view$after[1], view$after[2],
             col = colour, lwd = 2.5, length = 0.1)
      points(view$before[1], view$before[2], pch = 1, cex = 1.5, lwd = 2)
    }
    points(view$after[1], view$after[2], pch = 21, bg = colour, col = "black", cex = 1.5)
  }, res = 110, alt = "Animated Gibbs trajectory over joint posterior contours. Each move changes only one coordinate.")

  output$conditional_plot <- renderPlot({
    view <- current()
    distribution <- view$distribution
    limits <- range(distribution$quantile(c(0.001, 0.999)),
                    view$before[view$coordinate], view$after[view$coordinate])
    x <- sort(unique(c(seq(limits[1], limits[2], length.out = 400),
                       distribution$quantile(seq(0.001, 0.999, length.out = 400)))))
    y <- distribution$density(x)
    colour <- colours[view$coordinate]
    par(mar = c(4, 4, 2, 1))
    plot(x, y, type = "n", xlab = chain()$target$axes[view$coordinate],
         ylab = "Conditional density", main = paste("Update parameter", view$coordinate),
         ylim = c(0, max(y)))
    polygon(c(x[1], x, tail(x, 1)), c(0, y, 0),
            col = adjustcolor(colour, alpha.f = 0.15), border = NA)
    lines(x, y, col = colour, lwd = 2)
    abline(v = view$before[view$coordinate], col = "#666666", lty = 2)
    if (view$has_draw) abline(v = view$after[view$coordinate], col = colour, lwd = 2)
    legend("topright", c("Previous value", if (view$has_draw) "New draw"),
           col = c("#666666", if (view$has_draw) colour),
           lty = c(2, if (view$has_draw) 1), bty = "n", cex = 0.8)
  }, res = 110, alt = "Exact full conditional used by the current Gibbs update, with previous value and new draw marked.")

  output$step_table <- renderTable({
    view <- current()
    data.frame(Parameter = chain()$target$axes, Before = view$before,
               After = view$after,
               Status = ifelse(1:2 == view$coordinate,
                               if (view$has_draw) "Sampled" else "Next draw", "Held fixed"))
  }, digits = 4, rownames = FALSE)
  output$sample_count <- renderUI({
    p(sprintf("%d completed sweeps; burn-in setting: %d; retained so far: %d.",
              nrow(completed()), input$burn, nrow(retained())))
  })
  output$diagnostic_plot <- renderPlot({
    validate(need(nrow(retained()) >= 2, "Play beyond burn-in (at least two retained sweeps), or reduce burn-in."))
    samples <- retained()
    model <- chain()$target
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (j in 1:2) {
      trace <- c(chain()$states[1, j], completed()[, j])
      plot(0:nrow(completed()), trace, type = "l", col = colours[j],
           xlab = "Completed sweep", ylab = model$axes[j], main = paste("Trace:", model$axes[j]))
      abline(v = input$burn, lty = 2, col = "#666666")
      if (!is.null(model$mean)) abline(h = model$mean[j], lty = 3)
    }
    for (j in 1:2) {
      histogram <- hist(samples[, j], breaks = "FD", plot = FALSE)
      x <- seq(min(histogram$breaks), max(histogram$breaks), length.out = 400)
      reference <- if (!is.null(model$marginal)) model$marginal(x, j) else 0
      plot(histogram, freq = FALSE, col = "#DDDDDD", border = "white",
           xlab = model$axes[j], main = paste("Retained:", model$axes[j]),
           ylim = c(0, max(histogram$density, reference)))
      if (!is.null(model$marginal)) lines(x, reference, col = colours[j], lwd = 2)
    }
  }, res = 110, alt = "Trace plots and retained-sample marginal histograms for both parameters.")
  output$sample_summary <- renderTable({
    validate(need(nrow(retained()) >= 2, "At least two retained sweeps are needed for summaries."))
    result <- data.frame(Parameter = chain()$target$axes,
                         `Sample mean` = colMeans(retained()),
                         `Sample SD` = apply(retained(), 2, sd), check.names = FALSE)
    if (!is.null(chain()$target$mean)) {
      result[["Target mean"]] <- chain()$target$mean
      result[["Target SD"]] <- chain()$target$sd
    }
    result
  }, digits = 4, rownames = FALSE)
  output$acf_plot <- renderPlot({
    validate(need(nrow(retained()) >= 3, "At least three retained sweeps are needed for autocorrelation."))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    for (j in 1:2) {
      acf(retained()[, j], lag.max = min(40, nrow(retained()) - 1),
          main = paste("Autocorrelation:", chain()$target$axes[j]))
    }
  }, res = 110, alt = "Autocorrelation of each coordinate across retained complete sweeps.")

  output$model_equations <- renderUI({
    if (input$family == "poisson") {
      return(tagList(
        p("A genuine hierarchical posterior: n independent unit-exposure counts with total S. All Gamma parameters use rates, not scales."),
        tags$pre("x_i | lambda ~ Poisson(lambda)\nlambda | beta ~ Gamma(shape = a, rate = beta)\nbeta ~ Gamma(shape = c, rate = d)\n\nlambda | beta, x ~ Gamma(shape = a + S, rate = n + beta)\nbeta | lambda, x ~ Gamma(shape = a + c, rate = d + lambda)")
      ))
    }
    tagList(
      p("A single vector observation y has unknown location theta and known scale matrix S. ",
        "An improper flat location prior p(theta) proportional to 1 gives the proper joint posterior shown here. ",
        "For a Gaussian S is the covariance matrix; for Student-t its covariance is nu/(nu - 2) times S."),
      tags$pre("S = [[s1^2, rho*s1*s2], [rho*s1*s2, s2^2]]\ny | theta ~ Gaussian_2(theta, S) or Student-t_2(nu, theta, S)"),
      if (input$family == "gaussian") {
        tags$pre("theta1 | theta2, y ~ Normal(y1 + rho*s1/s2*(theta2-y2),\n                            variance = s1^2*(1-rho^2))\ntheta2 | theta1, y ~ Normal(y2 + rho*s2/s1*(theta1-y1),\n                            variance = s2^2*(1-rho^2))")
      } else {
        tags$pre("For coordinate j and the other coordinate k:\nm_j = y_j + rho*s_j/s_k*(theta_k-y_k)\ndelta_k = ((theta_k-y_k)/s_k)^2\n\ntheta_j | theta_k, y ~ Student-t(df = nu+1, location = m_j,\n  scale = s_j * sqrt((1-rho^2)*(nu+delta_k)/(nu+1)))")
      }
    )
  })
}

shinyApp(ui = ui, server = server)
