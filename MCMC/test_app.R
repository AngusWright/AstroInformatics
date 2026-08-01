# Run with: Rscript test_app.R
#
# This runner uses small chains so it can exercise every categorical UI
# combination quickly. It writes a detailed report to app-test-report.txt.

required_packages <- c(
  "shiny", "bslib", "mvtnorm", "bsicons", "magicaxis", "seqinr", "htmltools"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this test."
  )
}

source("app.R", local = TRUE)

report_path <- "app-test-report.txt"
test_cases <- expand.grid(
  model_type = c("gaussian", "binomial", "beta", "gamma"),
  prior_type = c("uniform", "gaussian"),
  sampler = c("gibbs", "mh"),
  plot_type = c("posterior", "parameters", "chain", "correlation"),
  pstyle = c(FALSE, TRUE),
  prior_vol = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

run_case <- function(case) {
  result <- tryCatch({
    shiny::testServer(server, {
      session$setInputs(
        model_type = case$model_type,
        num_data = 10,
        theta1 = 0.25,
        theta2 = 0.75,
        stdev = 0.5,
        rho = 0.5,
        prior = c(-2.5, 2.5),
        prior_type = case$prior_type,
        sampler = case$sampler,
        prop_sd = 0.2,
        n_chain = 1,
        n_samp = 2,
        n_burn = 1,
        step_thin = 1,
        plot_type = case$plot_type,
        pstyle = case$pstyle,
        prior_vol = case$prior_vol,
        annot_cex = 1.5
      )
      session$flushReact()
      stopifnot(
        is.list(output$plot),
        is.list(output$plot2),
        is.list(output$plot3),
        is.list(output$plot4)
      )
    })
    "PASS"
  }, error = function(error) {
    paste0("FAIL: ", conditionMessage(error))
  })

  sprintf(
    "%-4s model=%-8s prior=%-8s sampler=%-5s plot=%-11s lines=%-5s prior_volume=%-5s",
    if (identical(result, "PASS")) "PASS" else "FAIL",
    case$model_type,
    case$prior_type,
    case$sampler,
    case$plot_type,
    case$pstyle,
    case$prior_vol
  ) |> paste(result, sep = " | ")
}

results <- vapply(
  seq_len(nrow(test_cases)),
  function(index) run_case(test_cases[index, ]),
  character(1)
)
test_cases$passed <- startsWith(results, "PASS")
passed <- sum(test_cases$passed)
failed <- length(results) - passed
configuration_summary <- aggregate(
  passed ~ model_type + prior_type + sampler,
  data = test_cases,
  FUN = sum
)
configuration_totals <- aggregate(
  passed ~ model_type + prior_type + sampler,
  data = test_cases,
  FUN = length
)
names(configuration_totals)[4] <- "total"
configuration_summary <- merge(
  configuration_summary,
  configuration_totals,
  by = c("model_type", "prior_type", "sampler")
)
configuration_lines <- sprintf(
  "model=%-8s prior=%-8s sampler=%-5s %2d/%2d cases passed",
  configuration_summary$model_type,
  configuration_summary$prior_type,
  configuration_summary$sampler,
  configuration_summary$passed,
  configuration_summary$total
)

report <- c(
  "Shiny MCMC app test report",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Total cases: %d", length(results)),
  sprintf("Passed: %d", passed),
  sprintf("Failed: %d", failed),
  "",
  "Each case renders all four diagnostic outputs using one chain of 100 draws.",
  "",
  "Configuration summary:",
  configuration_lines,
  "",
  "Detailed results:",
  results
)

writeLines(report, report_path)
cat(paste(report[1:6], collapse = "\n"), "\n")
cat(sprintf("Detailed results: %s\n", normalizePath(report_path)))

if (failed > 0) {
  quit(status = 1)
}
