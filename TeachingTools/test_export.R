source("apps.R")
gallery <- paste(readLines("docs/index.html"), collapse = "\n")
overview <- paste(readLines("AstroInformatics.qmd"), collapse = "\n")
positions <- integer(length(apps))

verify_standalone <- function(payload) {
  directory <- tempfile("astro-export-check-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  for (name in c("app.R", "model.R", "teaching.R")) {
    writeLines(payload$content[match(name, payload$name)], file.path(directory, name))
  }
  app <- source(file.path(directory, "app.R"), local = new.env(), chdir = TRUE)$value
  stopifnot(inherits(app, "shiny.appobj"))
}

for (i in seq_along(apps)) {
  app <- apps[[i]]
  directory <- file.path("docs", app$path)
  stopifnot(file.exists(file.path(directory, "index.html")))
  link <- paste0('href="', app$path, '/"')
  positions[i] <- regexpr(link, gallery, fixed = TRUE)[1]
  stopifnot(positions[i] > 0,
            grepl(app$lecture, gallery, fixed = TRUE),
            grepl(paste0("](", app$path, "/)"), overview, fixed = TRUE))
  page <- paste(readLines(file.path(directory, "index.html")), collapse = "\n")
  stopifnot(grepl(paste0("<title>", htmltools::htmlEscape(app$title), "</title>"), page, fixed = TRUE))
  if (!is.null(app$static_files)) {
    asset_tags <- regmatches(page, gregexpr('(?:src|href)="[^"]+\\.(?:js|css)"', page, perl = TRUE))[[1]]
    references <- sub('^[^=]+="([^"]+)"$', "\\1", asset_tags)
    stopifnot(grepl(app$lecture, page, fixed = TRUE),
              !grepl("<script[^>]+(shinylive|webr)|app[.]json|https?://[^\" ]+[.](js|css)", page, ignore.case = TRUE),
              length(references) > 0, all(references %in% app$static_files),
              all(file.exists(file.path(directory, references))),
              !any(file.exists(file.path(directory, c("app.json", "edit/index.html")))))
    for (asset in app$static_files) {
      expected <- readBin(file.path(app$path, asset), "raw", n = file.info(file.path(app$path, asset))$size)
      actual <- readBin(file.path(directory, asset), "raw", n = file.info(file.path(directory, asset))$size)
      stopifnot(identical(actual, expected))
    }
    next
  }
  stopifnot(file.exists(file.path(directory, "app.json")), file.exists(file.path(directory, "edit/index.html")))
  payload <- jsonlite::fromJSON(file.path(directory, "app.json"))
  stopifnot(all(c("name", "content", "type") %in% names(payload)),
            "app.R" %in% payload$name,
            !any(grepl("Rhistory|DS_Store", payload$name)))
  if (app$shared) stopifnot(all(c("teaching.R", "model.R") %in% payload$name))
  for (j in seq_len(nrow(payload))) {
    if (payload$type[j] != "text") next
    source <- if (app$shared && payload$name[j] == "teaching.R") {
      "TeachingTools/teaching.R"
    } else file.path(app$path, payload$name[j])
    stopifnot(file.exists(source))
    expected <- sub("\n+$", "", paste(readLines(source, warn = FALSE), collapse = "\n"))
    actual <- sub("\n+$", "", gsub("\r\n", "\n", payload$content[j], fixed = TRUE))
    if (!identical(actual, expected)) stop(app$path, "/", payload$name[j], ": exported source differs")
  }
  if (app$shared) verify_standalone(payload)
}
stopifnot(all(diff(positions) > 0))
runtime <- "docs/shinylive/webr"
metadata <- readRDS(file.path(runtime, "packages/metadata.rds"))
for (package in metadata) {
  archive <- file.path(runtime, package$path)
  stopifnot(file.exists(archive), any(grepl("(^|/)DESCRIPTION$", utils::untar(archive, list = TRUE))))
}
cat(length(apps), "exported pages: catalogue order, lecture links, titles, source fidelity,\n",
    "standalone shared-module loading, runtime packages, source-view pages and local-history exclusion passed.\n")
