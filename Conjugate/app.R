library(shiny)
library(bslib)
source("conjugate.R", local = TRUE)

ui <- page_sidebar(
  title = "Conjugate priors: the role of data",
  fillable = FALSE,
  sidebar = sidebar(
    open = "open", width = 330,
    p(a("Read the lecture: Conjugate Priors and the Role of Data",
        href = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section3b_html.html#conjugate-priors")),
    p(a("All demonstrations", href = "../")),
    selectInput(
      "model", "Likelihood / conjugate prior",
      c("Normal / Normal (known variance)" = "normal",
        "Binomial / Beta" = "binomial",
        "Poisson / Gamma" = "poisson",
        "Gamma / Gamma (known shape)" = "gamma")
    ),
    checkboxInput("no_data", "Start with no data", FALSE),
    sliderInput("log_n", "Sample size n = round(10^x)",
                min = 0, max = 4, value = 1, step = 0.1),
    conditionalPanel(
      "input.model == 'normal'",
      sliderInput("sample_mean", "Observed sample mean",
                  min = -5, max = 5, value = 1, step = 0.1),
      sliderInput("sigma", "Known observation standard deviation (sigma)",
                  min = 0.1, max = 5, value = 2, step = 0.1),
      sliderInput("prior_mean", "Prior mean (b)",
                  min = -5, max = 5, value = -1, step = 0.1),
      sliderInput("prior_sd", "Prior standard deviation (s0)",
                  min = 0.05, max = 10, value = 1, step = 0.05)
    ),
    conditionalPanel(
      "input.model == 'binomial'",
      sliderInput("success_fraction", "Observed success fraction (rounded to s/n)",
                  min = 0, max = 1, value = 0.7, step = 0.01),
      sliderInput("beta_a", "Prior Beta shape (alpha)",
                  min = 0.5, max = 50, value = 2, step = 0.5),
      sliderInput("beta_b", "Prior Beta shape (beta)",
                  min = 0.5, max = 50, value = 5, step = 0.5)
    ),
    conditionalPanel(
      "input.model == 'poisson'",
      sliderInput("count_mean", "Observed counts per observation (rounded total)",
                  min = 0, max = 10, value = 4, step = 0.1)
    ),
    conditionalPanel(
      "input.model == 'gamma'",
      sliderInput("gamma_mean", "Observed sample mean",
                  min = 0.1, max = 10, value = 2, step = 0.1),
      sliderInput("shape", "Known observation shape (k)",
                  min = 0.5, max = 10, value = 3, step = 0.5)
    ),
    conditionalPanel(
      "input.model == 'poisson' || input.model == 'gamma'",
      sliderInput("gamma_a", "Prior Gamma shape (alpha)",
                  min = 0.5, max = 50, value = 3, step = 0.5),
      sliderInput("gamma_b", "Prior Gamma rate (beta, NOT scale)",
                  min = 0.1, max = 20, value = 2, step = 0.1)
    ),
    helpText("These controls specify observed sufficient statistics, not a new random sample. Changing the prior never changes the data."),
    hr(),
    h5("Gaussian experiments"),
    actionButton("strong_prior", "Prior dominates"),
    actionButton("many_data", "Data dominates"),
    actionButton("weak_prior", "Weak prior"),
    actionButton("reset", "Reset")
  ),
  p("Bayes' theorem: posterior is proportional to prior times likelihood. ",
    "With a conjugate pair, the prior and posterior belong to the same distribution family, ",
    "so the update is exact: no MCMC or density estimation is needed."),
  tabsetPanel(
    tabPanel(
      "Conjugate update",
      uiOutput("model_description"),
      plotOutput("density_plot", height = "390px"),
      p("The likelihood is scaled to unit area over the parameter for comparison. ",
        "It is not a probability distribution for the parameter unless a suitable flat prior is assumed. ",
        "Curves show central quantile ranges; divergent boundary densities are not drawn at the endpoints."),
      tableOutput("summary"),
      uiOutput("update_rule"),
      uiOutput("information")
    ),
    tabPanel(
      "The role of data",
      p("Hold the prior, observation noise (or known shape), and observed per-observation summary fixed, ",
        "then increase n. Integer count totals are rounded. This is a controlled comparison, ",
        "not a simulated sequence of new observations."),
      plotOutput("role_plot", height = "510px"),
      p("The shaded band is a central 95% credible interval for the unknown parameter, ",
        "not a prediction interval for a future observation. The vertical line marks the selected n (if nonzero)."),
      uiOutput("role_explanation")
    ),
    tabPanel(
      "What to try",
      h4("Prior knowledge versus data"),
      tags$ol(
        tags$li("Click 'Prior dominates': a small prior standard deviation means high prior precision. ",
                "With little data, the posterior stays near the prior, even when the data disagree."),
        tags$li("Click 'Data dominates': increasing n adds information, pulls the posterior towards ",
                "the sample estimate and reduces the influence of this fixed prior."),
        tags$li("Click 'Weak prior': a broad Normal prior has low precision. ",
                "The posterior approaches Normal(sample mean, sigma^2/n), even without a very large n."),
        tags$li("Select 'Start with no data': posterior equals prior. Then add observations and compare."),
        tags$li("Try the other pairs: conjugacy preserves the family, not the parameters. ",
                "Counts, exposure, or the sum of measurements update the prior analytically.")
      ),
      h4("Sufficiency: why only summaries are needed"),
      p("For Normal observations with known variance, n and the sample mean contain all data information ",
        "about the unknown mean. The individual measurements change the likelihood only by a ",
        "parameter-independent factor. For example, these two datasets give the same posterior ",
        "using the current Normal prior and known observation standard deviation:"),
      tableOutput("sufficiency"),
      h4("Important qualifications"),
      p("For the Normal model, posterior precision is prior precision plus n times observation precision. ",
        "Consequently its variance always decreases as n increases. This exact precision-addition ",
        "statement is specific to the known-variance Normal example, not every conjugate family."),
      p("With a fixed proper prior, sufficiently informative data dominate in these examples. ",
        "Convergence to a true parameter also requires an appropriate sampling model and prior support ",
        "near that truth. Data cannot rescue a parameter region excluded by the prior."),
      p("A Normal prior standard deviation tending to infinity is a flat, improper-prior limit, ",
        "not a proper uniform probability distribution on the real line. All priors selectable here are proper."),
      p("Conjugacy is a computational convenience, not a reason to choose a prior that misrepresents your knowledge."),
      a("Lecture notes: Conjugate Priors and The Role of Data",
        href = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section3b_html.html#conjugate-priors")
    )
  )
)

server <- function(input, output, session) {
  sample_size <- reactive({
    if (isTRUE(input$no_data)) 0 else round(10^input$log_n)
  })

  analysis_at <- function(n, model = input$model) {
    statistic <- switch(
      model,
      normal = input$sample_mean,
      binomial = round(n * input$success_fraction),
      poisson = round(n * input$count_mean),
      gamma = n * input$gamma_mean
    )
    prior <- switch(
      model,
      normal = c(input$prior_mean, input$prior_sd),
      binomial = c(input$beta_a, input$beta_b),
      c(input$gamma_a, input$gamma_b)
    )
    conjugate_update(model, n, statistic, prior[1], prior[2],
                     sigma = input$sigma, shape = input$shape)
  }
  analysis <- reactive(analysis_at(sample_size()))

  normal_preset <- function(log_n, prior_sd) {
    updateSelectInput(session, "model", selected = "normal")
    updateCheckboxInput(session, "no_data", value = FALSE)
    updateSliderInput(session, "log_n", value = log_n)
    updateSliderInput(session, "sample_mean", value = 1)
    updateSliderInput(session, "sigma", value = 2)
    updateSliderInput(session, "prior_mean", value = -1)
    updateSliderInput(session, "prior_sd", value = prior_sd)
  }
  observeEvent(input$strong_prior, normal_preset(0, 0.2))
  observeEvent(input$many_data, normal_preset(4, 0.2))
  observeEvent(input$weak_prior, normal_preset(1, 10))
  observeEvent(input$reset, normal_preset(1, 1))

  output$model_description <- renderUI({
    descriptions <- c(
      normal = "Infer the mean theta of Normal observations with known variance sigma^2. The conjugate prior and posterior are Normal.",
      binomial = "Infer the success probability theta from s successes in n independent trials. The conjugate prior and posterior are Beta.",
      poisson = "Infer a Poisson rate theta from n independent, equal unit-exposure counts. The conjugate prior and posterior are Gamma.",
      gamma = "Infer the rate theta of Gamma observations with known shape k. The conjugate prior and posterior are Gamma."
    )
    tagList(h4(descriptions[[input$model]]), p(analysis()$data_description))
  })

  output$density_plot <- renderPlot({
    result <- analysis()
    distributions <- list(result$prior, result$likelihood, result$posterior)
    colours <- c("#0072B2", "#666666", "#D55E00")
    labels <- c("Prior", "Likelihood (unit area)", "Posterior")
    active <- which(!vapply(distributions, is.null, logical(1)))
    grids <- lapply(distributions[active], function(d) {
      d$quantile(seq(0.001, 0.999, length.out = 1200))
    })
    limits <- range(unlist(grids))
    x <- sort(unique(c(seq(limits[1], limits[2], length.out = 1200),
                       unlist(grids))))
    # Do not evaluate Beta/Gamma densities at singular support boundaries.
    if (input$model == "binomial") x <- x[x > 0 & x < 1]
    if (input$model %in% c("poisson", "gamma")) x <- x[x > 0]
    curves <- lapply(distributions[active], function(d) d$density(x))
    par(mar = c(4, 4, 1, 1))
    plot(range(x), c(0, max(unlist(curves))), type = "n",
         xlab = if (input$model %in% c("poisson", "gamma")) "Unknown rate theta" else "Unknown parameter theta",
         ylab = "Density")
    for (j in seq_along(active)) {
      lines(x, curves[[j]], col = colours[active[j]],
            lty = active[j], lwd = 2.5)
    }
    legend("topright", legend = labels[active], col = colours[active],
           lty = active, lwd = 2.5, bty = "n", cex = 0.85)
  }, res = 110,
  alt = "Prior, area-normalised likelihood, and posterior densities over the unknown parameter.")

  output$summary <- renderTable({
    result <- analysis()
    distributions <- list(Prior = result$prior, Posterior = result$posterior)
    intervals <- lapply(distributions, function(d) d$quantile(c(0.025, 0.975)))
    data.frame(
      Distribution = names(distributions),
      Family = vapply(distributions, `[[`, character(1), "label"),
      Mean = vapply(distributions, `[[`, numeric(1), "mean"),
      `Standard deviation` = vapply(distributions, `[[`, numeric(1), "sd"),
      `95% credible interval` = vapply(intervals, function(q) {
        sprintf("[%.4g, %.4g]", q[1], q[2])
      }, character(1)),
      check.names = FALSE
    )
  }, digits = 4, striped = TRUE, rownames = FALSE)

  output$update_rule <- renderUI({
    rules <- list(
      normal = tagList(
        p("Write d = 1/s0^2 and tau = 1/sigma^2 (precisions)."),
        tags$pre("Prior:     theta ~ Normal(b, variance = 1/d)\nPosterior: theta | x ~ Normal((d*b + n*tau*xbar)/(d + n*tau),\n                              variance = 1/(d + n*tau))")
      ),
      binomial = tags$pre("Prior:     theta ~ Beta(alpha, beta)\nData:      s | theta ~ Binomial(n, theta)\nPosterior: theta | s ~ Beta(alpha + s, beta + n - s)"),
      poisson = tags$pre("Prior:     theta ~ Gamma(shape = alpha, rate = beta)\nData:      x_i | theta ~ Poisson(theta)\nPosterior: theta | x ~ Gamma(shape = alpha + sum(x), rate = beta + n)"),
      gamma = tags$pre("Prior:     theta ~ Gamma(shape = alpha, rate = beta)\nData:      x_i | theta ~ Gamma(shape = k, rate = theta)\nPosterior: theta | x ~ Gamma(shape = alpha + n*k, rate = beta + sum(x))")
    )
    tagList(h4("Exact conjugate update"), rules[[input$model]],
            p("Normal formulas use variances; Gamma formulas always use rates, not scales."))
  })

  output$information <- renderUI({
    result <- analysis()
    if (result$n == 0) {
      return(p(strong("No observations: "), "the posterior is exactly the prior; there is no data estimate or likelihood curve."))
    }
    tagList(
      p(sprintf("Posterior mean = %.1f%% prior mean + %.1f%% data estimate (%.4g).",
                100 * result$prior_weight, 100 * (1 - result$prior_weight),
                result$estimate)),
      if (input$model == "normal") {
        p(sprintf("Posterior precision = %.4g + %.4g = %.4g. Posterior standard deviation = %.4g.",
                  1 / input$prior_sd^2, result$n / input$sigma^2,
                  1 / result$posterior$sd^2, result$posterior$sd))
      },
      p(result$weight_description)
    )
  })

  output$role_plot <- renderPlot({
    ns <- sort(unique(c(1, round(10^seq(0, 4, length.out = 160)))))
    results <- lapply(ns, analysis_at)
    means <- vapply(results, function(r) r$posterior$mean, numeric(1))
    intervals <- vapply(results, function(r) r$posterior$quantile(c(0.025, 0.975)),
                        numeric(2))
    weights <- vapply(results, `[[`, numeric(1), "prior_weight")
    estimates <- vapply(results, `[[`, numeric(1), "estimate")
    width <- session$clientData$output_role_plot_width
    if (is.null(width) || width >= 700) par(mfrow = c(1, 2)) else par(mfrow = c(2, 1))
    par(mar = c(4, 4, 2, 1))
    plot(ns, means, type = "n", log = "x", ylim = range(intervals, estimates, analysis()$prior$mean),
         xlab = "Number of observations n", ylab = "Unknown parameter theta",
         main = "Posterior mean and uncertainty")
    polygon(c(ns, rev(ns)), c(intervals[1, ], rev(intervals[2, ])),
            col = adjustcolor("#D55E00", alpha.f = 0.18), border = NA)
    abline(h = analysis()$prior$mean, col = "#0072B2", lty = 2, lwd = 2)
    lines(ns, estimates, col = "#666666", lty = 3, lwd = 2)
    lines(ns, means, col = "#D55E00", lwd = 2)
    if (sample_size() > 0) abline(v = sample_size(), lty = 3)
    legend("topright", c("Posterior mean", "Prior mean", "Data estimate"),
           col = c("#D55E00", "#0072B2", "#666666"), lty = 1:3,
           lwd = 2, bty = "n", cex = 0.8)
    plot(ns, weights, type = "l", log = "x", ylim = c(0, 1), lwd = 2,
         col = "#0072B2", xlab = "Number of observations n",
         ylab = "Weight of prior mean", main = "Relative prior influence")
    if (sample_size() > 0) abline(v = sample_size(), lty = 3)
  }, res = 110,
  alt = "Posterior mean and 95 percent credible interval, and prior weight, as sample size increases.")

  output$role_explanation <- renderUI({
    tagList(
      p(analysis()$weight_description),
      if (input$model == "normal") {
        p("The posterior mean is gamma_n*b + (1 - gamma_n)*xbar, where ",
          "gamma_n = d/(d + n*tau). More data or a smaller known sigma increase data precision. ",
          "A larger prior standard deviation decreases prior precision. ",
          "For fixed prior and sigma, as n grows the posterior approaches Normal(xbar, sigma^2/n).")
      } else {
        p("Here the weight is the exact coefficient in a weighted-mean decomposition, ",
          "not a universal fraction of information or probability that the prior is correct. ",
          "For Gamma rate inference, if the observed sample mean is extremely small, ",
          "the prior rate can remain influential even at large n.")
      }
    )
  })

  output$sufficiency <- renderTable({
    centre <- input$sample_mean
    datasets <- list(centre + c(-1, -0.5, 0.5, 1),
                     centre + c(-2, -1, 1, 2))
    results <- lapply(datasets, function(x) {
      conjugate_update("normal", length(x), mean(x), input$prior_mean,
                       input$prior_sd, sigma = input$sigma)
    })
    data.frame(
      Dataset = vapply(datasets, function(x) paste(format(x, trim = TRUE), collapse = ", "), character(1)),
      n = 4, `Sample mean` = centre,
      `Posterior mean` = vapply(results, function(r) r$posterior$mean, numeric(1)),
      `Posterior SD` = vapply(results, function(r) r$posterior$sd, numeric(1)),
      check.names = FALSE
    )
  }, digits = 4, rownames = FALSE)
}

shinyApp(ui = ui, server = server)
