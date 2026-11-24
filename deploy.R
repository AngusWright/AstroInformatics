source("apps.R")
output_dir <- "docs"
escape <- function(text) as.character(htmltools::htmlEscape(text))

export_lesson <- function(app) {
  directory <- app$path
  if (!is.null(app$static_files)) {
    destination <- file.path(output_dir, app$path)
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    message("Copying static page ", app$path)
    for (asset in app$static_files) {
      target <- file.path(destination, asset)
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      stopifnot(file.copy(file.path(directory, asset), target, overwrite = TRUE))
    }
    stale <- file.path(destination, c("app.json", "edit/index.html"))
    unlink(stale)
    stopifnot(!any(file.exists(stale)))
    return(invisible(NULL))
  }
  if (app$shared) {
    staging <- tempfile(paste0("astro-", app$path, "-"))
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE))
    sources <- list.files(directory, full.names = TRUE)
    stopifnot(all(file.copy(sources, staging, recursive = TRUE)),
              file.copy("TeachingTools/teaching.R", file.path(staging, "teaching.R")))
    directory <- staging
  }
  message("Exporting ", app$path)
  shinylive::export(directory, output_dir, subdir = app$path)
  for (relative in c("index.html", "edit/index.html")) {
    page <- file.path(output_dir, app$path, relative)
    markup <- readLines(page)
    if (relative == "index.html") {
      markup <- sub("<title>[^<]*</title>", paste0("<title>", escape(app$title), "</title>"), markup)
    }
    writeLines(sub("[[:blank:]]+$", "", markup), page)
  }
}
invisible(lapply(apps, export_lesson))

groups <- unique(vapply(apps, `[[`, character(1), "group"))
cards <- vapply(groups, function(group) {
  selected <- apps[vapply(apps, function(app) app$group == group, logical(1))]
  entries <- vapply(seq_along(selected), function(i) {
    app <- selected[[i]]
    paste0(
      '<article class="card"><span class="eyebrow">',
      if (group == "Priority demonstrations") sprintf("START HERE / %02d", i) else escape(group),
      '</span><h3><a href="', app$path, '/">', escape(app$title), '</a></h3>',
      '<p>', escape(app$description), '</p><div class="card-links"><a class="explore" href="',
      app$path, '/">Explore demonstration &rarr;</a><a href="', escape(app$lecture),
      '">Lecture: ', escape(app$section), '</a></div></article>'
    )
  }, character(1))
  paste0('<section><h2>', escape(group), '</h2><div class="grid">',
         paste(entries, collapse = "\n"), '</div></section>')
}, character(1))

writeLines(
  c(
    "<!doctype html>",
    '<html lang="en">',
    "  <head>",
    '    <meta charset="utf-8">',
    '    <meta name="viewport" content="width=device-width, initial-scale=1">',
    "    <title>Astroinformatics | See the statistics</title>",
    '    <meta name="description" content="Interactive astronomical statistics: selection effects, Bayesian inference, sampling, simulation, and statistical foundations.">',
    "    <style>",
    "      *{box-sizing:border-box}body{margin:0;background:#f5f7fb;color:#25344b;font:16px/1.6 system-ui,-apple-system,sans-serif}",
    "      header{background:linear-gradient(125deg,#182a46,#255e73);color:white;padding:64px max(24px,calc((100vw - 1180px)/2))}",
    "      header p{max-width:760px;color:#dce8f1;font-size:1.15rem}h1{font-size:clamp(2rem,5vw,3.8rem);line-height:1.1;margin:12px 0}",
    "      header a{color:#c6eee9}.eyebrow{font-size:.73rem;letter-spacing:.1em;text-transform:uppercase;font-weight:700;color:#087f8c}",
    "      header .eyebrow{color:#b1e6df}main{max-width:1228px;margin:auto;padding:16px 24px 60px}section{margin:36px 0}",
    "      h2{font-size:1.45rem;margin:0 0 18px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,310px),1fr));gap:20px}",
    "      .card{background:white;border:1px solid #e0e5ef;border-radius:16px;padding:24px;display:flex;flex-direction:column;box-shadow:0 5px 16px #25344b05}",
    "      .card h3{line-height:1.25;font-size:1.25rem;margin:12px 0}.card h3 a{color:#25344b;text-decoration:none}.card p{color:#536176;margin:0 0 24px}",
    "      a{color:#365b99;text-underline-offset:3px}a:hover{text-decoration:underline}.card-links{display:grid;gap:12px;margin-top:auto;font-size:.84rem}",
    "      .explore{font-size:.95rem;font-weight:700}footer{border-top:1px solid #dfe5ef;padding:24px;color:#536176;text-align:center}",
    "      a:focus-visible{outline:3px solid #d56534;outline-offset:4px}details{max-width:780px}summary{cursor:pointer;font-weight:600}",
    "    </style>",
    "  </head>",
    "  <body>",
    '    <header><span class="eyebrow">Astroinformatics / Interactive lecture companions</span>',
    '      <h1>See the statistics.</h1>',
    '      <p>Explore the principles behind astronomical inference, one clear visual experiment at a time. Change a few assumptions, reveal the hidden truth, and follow the evidence.</p>',
    '      <a href="https://anguswright.github.io/LectureNotes/IntroductionToStatistics_materials.html">Open the lecture series</a></header>',
    "    <main>",
    '      <p>Start with the numbered priority demonstrations, or revisit a familiar inference tool. Every page includes its relevant HTML lecture notes, a guided experiment, and the assumptions behind the model.</p>',
    cards,
    '      <details><summary>About these demonstrations</summary><p>These are controlled teaching simulations, not analyses of real survey data. Mathematical details and caveats are kept in each page&rsquo;s explanation view. The galaxy-encounter page is a restricted statistical companion to the original Holmberg simulation prototype, not a self-consistent N-body calculation.</p></details>',
    "    </main>",
    '    <footer>Astroinformatics &middot; Statistical principles for astronomers and physicists</footer>',
    "  </body>",
    "</html>"
  ),
  file.path(output_dir, "index.html")
)
