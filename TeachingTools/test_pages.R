root <- if (dir.exists("TeachingTools")) "." else ".."
arguments <- commandArgs(trailingOnly = TRUE)
directories <- list.dirs(root, recursive = FALSE, full.names = TRUE)
pages <- directories[vapply(directories, function(directory) {
  file <- file.path(directory, "app.R")
  file.exists(file) && any(grepl("lesson_components(spec", readLines(file, warn = FALSE), fixed = TRUE))
}, logical(1))]
catalogue <- new.env()
sys.source(file.path(root, "apps.R"), envir = catalogue)
registered <- Filter(function(app) app$shared, catalogue$apps)
stopifnot(setequal(basename(pages), names(registered)))
if (length(arguments)) pages <- pages[basename(pages) %in% arguments]
stopifnot(length(pages) > 0)

for (directory in pages) {
  app <- new.env()
  source(file.path(directory, "app.R"), local = app, chdir = TRUE)
  config <- app$spec
  defaults <- setNames(lapply(app$spec$controls, `[[`, "value"),
                       vapply(app$spec$controls, `[[`, character(1), "id"))
  defaults <- c(defaults, list(seed = 42, delay = "200"))
  stopifnot(grepl("^https://anguswright.github.io/LectureNotes/", app$spec$lecture))
  stopifnot(identical(app$spec$lecture, registered[[basename(directory)]]$lecture))
  html <- htmltools::renderTags(app$ui)$html
  stopifnot(grepl(app$spec$lecture, html, fixed = TRUE),
            grepl("The principle:", html, fixed = TRUE))
  shiny::testServer(app$server, {
    do.call(session$setInputs, defaults)
    original <- result()
    stopifnot(is.list(original), !playing(),
              nchar(output$plot$src) > 0, nchar(output$metrics$html) > 0)
    if (config$steps > 0) {
      session$setInputs(rewind = 1)
      stopifnot(frame() == 0, nchar(output$plot$src) > 0)
      session$setInputs(play = 1)
      stopifnot(playing(), frame() == 1)
      session$elapse(201)
      stopifnot(frame() == 2)
      session$setInputs(play = 2)
      paused <- frame()
      session$elapse(1000)
      stopifnot(!playing(), frame() == paused)
      session$setInputs(step = 1)
      stopifnot(frame() == paused + 1)
      session$setInputs(finish = 1)
      stopifnot(frame() == config$steps, !playing(),
                isTRUE(all.equal(result(), original)))
    }
    for (control in config$controls) {
      values <- if (control$type == "choice") unname(control$choices) else c(control$min, control$max)
      for (value in values) {
        do.call(session$setInputs, setNames(list(value), control$id))
        stopifnot(!playing(), nchar(output$plot$src) > 0,
                  nchar(output$metrics$html) > 0)
      }
      do.call(session$setInputs, setNames(list(control$value), control$id))
    }
    for (preset in config$presets) {
      do.call(session$setInputs, preset)
      stopifnot(!playing(), nchar(output$plot$src) > 0, nchar(output$metrics$html) > 0)
    }
    do.call(session$setInputs, defaults)
    stopifnot(isTRUE(all.equal(result(), original)))
    session$setInputs(redraw = 1)
    stopifnot(generation() == 1L, !playing(), nchar(output$plot$src) > 0)
  })
  cat(basename(directory), ": UI, models, control boundaries, reproducibility, plots and playback passed.\n")
}
