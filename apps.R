apps <- jsonlite::fromJSON("html-tools.json", simplifyVector = FALSE)
apps <- lapply(apps, function(app) {
  app$static_files <- unlist(app$static_files, use.names = FALSE)
  # This flag describes the retained R reference's use of teaching.R, not HTML.
  app$shared <- !app$path %in% c("Conjugate", "Gibbs", "MCMC", "ShiftFitting")
  app
})
names(apps) <- vapply(apps, `[[`, character(1), "path")
